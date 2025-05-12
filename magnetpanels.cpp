/*
*  magnetpanels.cpp
*  
*  Created by Duncan MacGregor
*  University of Edinburgh 2022
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*
*/

#include "magnetmodel.h"


MagNetBox::MagNetBox(MagNetModel *magnetmodel, MainFrame *main, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(magnetmodel, title, pos, size, "MagNet", 0, 1)
{
	column = 0;
	mod = magnetmodel;

	InitMenu();

	//SetModFlag(ID_Random, "randomflag", "Random Run", 1); 
	//SetModFlag(ID_seedgen, "netgen", "Net Gen", 1); 
	//SetModFlag(ID_singletrans, "singletrans", "Single MagNet", 0); 	
	//SetModFlag(ID_netdiag, "netdiag", "Net Diagnostic", 0); 

	SetModFlag(ID_spikemode, "spikemode", "Spiking Mod", 1);
	SetModFlag(ID_secmode, "secmode", "Secretion/Plasma Mod", 0);
	SetModFlag(ID_plasmamode, "plasmamode", "Plasma Mod", 0);
	SetModFlag(ID_inputgen, "inputgen", "Input Gen", 0); 
	SetModFlag(ID_realtime, "realtime", "Real Time", 0); 
	SetModFlag(ID_analysis, "netanalysis", "Net Analysis", 0); 


	// Parameter controls
	//
	// AddCon(tag string, display string, initial value, click increment, decimal places)
	// ----------------------------------------------------------------------------------

	paramset.AddCon("runtime", "Run Time", 2000, 1, 0);
	paramset.AddCon("numneurons", "Neurons", 10, 1, 0);
	paramset.AddCon("numruns", "Num Runs", 1, 1, 0);
	paramset.AddCon("netrate", "Net Rate", 100, 1, 0);  // the bigger the less accurate but faster. Minimum is 1 (synchronizing threads every ms)
	paramset.AddCon("osmorate", "Osmo Rate", 100, 1, 0); 
	paramset.AddCon("osmo_hstep", "Osmo hStep", 1, 1, 0); 
	paramset.AddCon("buffrate", "Buff Rate", 1000, 1, 0); 
	paramset.AddCon("synvarsd", "SynVar SD", 0, 0.05, 2);
	paramset.AddCon("inputcells", "inputcells", 200, 1, 0); 
	paramset.AddCon("neurosyn", "neurosyn", 100, 1, 0);
	paramset.AddCon("netinput", "Net Input", 100, 10, 2);
	paramset.AddCon("netIratio", "Net Iratio", 0.5, 0.1, 2);
	paramset.AddCon("popscale", "Pop Scale", 1000, 10, 2);
	paramset.AddCon("disprate", "Disp Rate", 1000, 10, 0);
	paramset.AddCon("storeinit", "Store Init", 2000000, 100000, 0);
	paramset.AddCon("mRNAinit", "mRNAinit", 20, 1, 2); 

	paramset.GetCon("runtime")->SetMinMax(0, 2000000);

	ParamLayout(2);

	SetPanel(ID_Grid, mod->gridbox);
	SetPanel(ID_Protocol, mod->protobox);
	SetPanel(ID_Signal, mod->signalbox);
	SetPanel(ID_Spike, mod->spikebox);
	SetPanel(ID_Net, mod->neurodatabox);

	// ----------------------------------------------------------------------------------

	defbutt = 1;
	wxBoxSizer *runbox = RunBox();

	wxBoxSizer *paramfilebox = StoreBoxSync();
	AddButton(ID_Compare, "Comp", 40, paramfilebox);

	wxBoxSizer *seedbox = new wxBoxSizer(wxHORIZONTAL);
	paramset.AddNum("modseed", "", 0, 0, 0, 80);
	paramset.GetCon("modseed")->SetMinMax(0, 1000000000000);
	seedbox->Add(paramset.GetCon("modseed"), 0, wxALIGN_CENTRE_VERTICAL);
	seedcheck = SetModCheck(ID_seedcheck, "seedgen", "Seed", true);
	seedbox->Add(seedcheck, 0, wxALIGN_CENTRE_VERTICAL|wxLEFT, 5);
	initcheck = SetModCheck(ID_initcheck, "netinit", "Init", true);
	seedbox->Add(initcheck, 0, wxALIGN_CENTRE_VERTICAL|wxLEFT, 5);
	storecheck = SetModCheck(ID_storereset, "storereset", "Res", true);
	seedbox->Add(storecheck, 0, wxALIGN_CENTRE_VERTICAL|wxLEFT, 5);


	buttonbox = new wxBoxSizer(wxHORIZONTAL);

	AddButton(ID_Spike, "Neuro", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();
	AddButton(ID_Signal, "Signal", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();
	AddButton(ID_Protocol, "Proto", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();
	AddButton(ID_Grid, "Grid", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();
	AddButton(ID_Net, "Net", buttonwidth, buttonbox);

	//synccheck = new wxCheckBox(panel, wxID_ANY, "Sync");
	//synccheck->SetValue(true);
	//paramfilebox->Add(synccheck, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);

	//buttonbox = new wxBoxSizer(wxHORIZONTAL);
	//paramfilebox->AddSpacer(5);
	//AddButton(ID_Output, "Output", 50, paramfilebox);

	status = StatusBar();
	wxBoxSizer *statusbox = new wxBoxSizer(wxHORIZONTAL);
	statusbox->Add(status, 1, wxEXPAND);

	mainbox->AddSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer();
	mainbox->Add(runbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	mainbox->AddSpacer(5);
	mainbox->Add(seedbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	//mainbox->AddSpacer(5);
	mainbox->Add(buttonbox, 0, wxALIGN_CENTRE_HORIZONTAL | wxALIGN_CENTRE_VERTICAL | wxALL, 0);
	mainbox->Add(paramfilebox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	mainbox->AddSpacer(2);
	//mainbox->AddStretchSpacer();
	mainbox->Add(statusbox, 0, wxEXPAND);
	//mainbox->AddStretchSpacer();
	//mainbox->Add(buttonbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);	

	panel->Layout();

	Connect(ID_paramstore, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetBox::OnParamStore));
	Connect(ID_paramload, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetBox::OnParamLoad));
	Connect(ID_Run, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetBox::OnRun));
	Connect(ID_Progress, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(MagNetBox::OnProgress));
	Connect(ID_Plot, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(MagNetBox::OnPlot));
	Connect(ID_Compare, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetBox::OnParamLoad));
}


void MagNetBox::OnPlot(wxCommandEvent& event)
{
	int count = event.GetInt();
	
	if((*modflags)["realtime"]) {
		mod->neurodatabox->NeuroData();
		mod->mainwin->scalebox->GraphUpdate();
	}
}


void MagNetBox::OnProgress(wxCommandEvent& event)
{
	int count = event.GetInt();
	SetCount(count);
}


void MagNetBox::OnParamStore(wxCommandEvent& event)
{
	wxString filetag, filepath, text;
	int i;

	if(synccheck->GetValue()) {
		filetag = paramstoretag->GetValue();
		mod->spikebox->StoreParam(filetag);
		mod->secbox->StoreParam(filetag);
		mod->signalbox->StoreParam(filetag);
		mod->protobox->StoreParam(filetag);
		mod->synthbox->StoreParam(filetag);
	}
	ParamBox::OnParamStore(event);	

	// Neuron Initialisation Store
	if(initcheck->GetValue() == false) {
		diagbox->Write("MagNetBox writing neuroinit file\n");
		filepath = mod->GetPath() + "/Init";
		if(!wxDirExists(filepath)) wxMkdir(filepath);
		filetag = paramstoretag->GetValue();
		TextFile neurofile;
		neurofile.New(filepath + "/" + filetag + "-neuroinit.txt");

		for(i=0; i<mod->numneurons; i++) {
			text.Printf("neuro %d  mRNAinit %.4f  store %.4f  synvar %.4f", i, mod->modneurons[i].mRNAinit, mod->modneurons[i].storeinit, mod->modneurons[i].synvar);
			neurofile.WriteLine(text);
		}
		neurofile.Close();
	}
}


void MagNetBox::OnParamLoad(wxCommandEvent& event)
{
	wxString filetag, filepath;
	bool compmode = false;
	bool check;
	wxString readline, tag, text;
	double numdat;
	int neurodex;

	if(event.GetId() == ID_Compare) compmode = true;

	if(synccheck->GetValue()) {
		filetag = paramstoretag->GetValue();
		mod->spikebox->ParamLoad(filetag, compmode);
		mod->secbox->ParamLoad(filetag, compmode);
		mod->signalbox->ParamLoad(filetag, compmode);
		mod->protobox->ParamLoad(filetag, compmode);
		mod->synthbox->ParamLoad(filetag, compmode);
	}
	ParamLoad(filetag, compmode);	

	// Neuron Initialisation Load
	if(initcheck->GetValue() == false) {
		diagbox->Write("MagNetBox loading neuroinit file\n");
		filetag = paramstoretag->GetValue();
		filepath = mod->GetPath() + "/Init";
		TextFile neurofile;
		check = neurofile.Open(filepath + "/" + filetag + "-neuroinit.txt");
		if(!check) mainwin->diagbox->Write("neuroinit file not found\n");
		else {
			readline = neurofile.ReadLine();
			while(!readline.IsEmpty()) {
				//diagbox->Write("neuroinit read " + readline + "\n");
				neurodex = ParseLong(&readline, 'o');
				mod->modneurons[neurodex].mRNAinit = ParseDouble(&readline, 't');
				mod->modneurons[neurodex].storeinit = ParseDouble(&readline, 'e');
				mod->modneurons[neurodex].synvar = ParseDouble(&readline, 'r');
				
				diagbox->Write(text.Format("MagNetBox neuron %d loaded %.4f %.2f %.4f\n", 
					neurodex, mod->modneurons[neurodex].mRNAinit, mod->modneurons[neurodex].storeinit, mod->modneurons[neurodex].synvar));
				readline = neurofile.ReadLine(); 

				//mod->modneurons[neurodex].initflag = true;
				//(*mod->modneurons[neurodex].spikeparams)["synvar"] = mod->modneurons[neurodex].synvar;
				//(*mod->modneurons[neurodex].synthparams)["mRNAinit"] = mod->modneurons[neurodex].mRNAinit;
				//(*mod->modneurons[neurodex].secparams)["Rinit"] = mod->modneurons[neurodex].storeinit;
			}
		}
		neurofile.Close();
	}

	// Set neuron count
	ParamStore *netparams = GetParams();
	mod->numneurons = (*netparams)["numneurons"];
	mod->NeuroGen();
	mod->neurodatabox->neurocount = mod->numneurons;
	mod->neurodatabox->NeuroData();

}


void MagNetBox::OnRun(wxCommandEvent& event)
{
	(*mod->modeflags)["prototype"] = 0;

	countmark = 0;
	mod->RunModel();
}


MagSynthBox::MagSynthBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(mod, title, pos, size, "MAGSYNTH")
{
	int labelwidth = 60;
	column = 0;
	int panelmode = 0;

	panelmode = 1;  // Parameters names for vaso synth paper

	InitMenu();

	SetModFlag(ID_boxflag, "boxflag", "Use Box", 0); 
	SetModFlag(ID_VSflag, "VSflag", "V Synth Store", 0); 
	SetModFlag(ID_transflag, "transflag", "Spike dependent translate", 0); 
	SetModFlag(ID_synthmode, "synthmode", "Synth Mode", 1); 
	SetModFlag(ID_transmode, "transmode", "Trans Mode", 1); 
	SetModFlag(ID_scalemode, "scalemode", "Synth Scale", 1); 
	SetModFlag(ID_polymode, "polymode", "Poly(A) mRNA", 1); 
	SetModFlag(ID_decaymode, "decaymode", "mRNA decay", 0); 

	if(panelmode == 1) {
		paramset.AddCon("mRNAinit", "m_init", 20, 1, 2, labelwidth); 
		paramset.AddCon("mRNAmax", "m_max", 100, 1, 2, labelwidth); 
		paramset.AddCon("mRNAhalflife", "m HL", 1000, 1, 1, labelwidth); 
		paramset.AddCon("synthdel", "synthdel", 0, 360, 0, labelwidth);
		paramset.AddCon("halflifeTS", "T HL", 1000, 1, 1, labelwidth); 
		paramset.AddCon("kTS", "kT", 1, 0.01, 4, labelwidth); 
		paramset.AddCon("maxTL", "maxTL", 5, 0.01, 4, labelwidth); 
		paramset.AddCon("kTL", "kTL", 1, 0.01, 6, labelwidth); 
		paramset.AddCon("halflifeTL", "TL HL", 1000, 1, 1, labelwidth); 
		paramset.AddCon("basalTL", "TL_basal", 1, 0.01, 4, labelwidth); 
		paramset.AddCon("rateSR", "s_r", 0.01, 0.001, 4, labelwidth); 
		paramset.AddCon("synscale", "s_scale", 0.0001, 0.00, 7, labelwidth);

		paramset.SetMinMax("halflifeVS", 0, 500000);
	}
	else {
		paramset.AddCon("mRNAinit", "mRNAinit", 20, 1, 2, labelwidth); 
		paramset.AddCon("mRNAmax", "mRNAmax", 100, 1, 2, labelwidth); 
		paramset.AddCon("mRNAhalflife", "mRNA HL", 1000, 1, 1, labelwidth); 
		paramset.AddCon("vtrans", "vtrans", 0.01, 0.001, 6, labelwidth);
		paramset.AddCon("synthdel", "synthdel", 0, 360, 0, labelwidth);
		paramset.AddCon("halflifeTS", "halflifeTS", 1000, 1, 1, labelwidth); 
		paramset.AddCon("kTS", "kTS", 1, 0.01, 4, labelwidth); 
		paramset.AddCon("maxTL", "maxTL", 5, 0.01, 4, labelwidth); 
		paramset.AddCon("kTL", "kTL", 1, 0.01, 6, labelwidth); 
		paramset.AddCon("halflifeTL", "halflifeTL", 1000, 1, 1, labelwidth); 
		paramset.AddCon("basalTL", "basalTL", 1, 0.01, 4, labelwidth); 
		paramset.AddCon("rateSR", "rateSR", 0.01, 0.001, 4, labelwidth); 
		paramset.AddCon("synscale", "synscale", 0.0001, 0.00, 7, labelwidth);

		paramset.SetMinMax("halflifeVS", 0, 500000);
	}

	

	ParamLayout(2);

	mainbox->AddStretchSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer(10);
	panel->Layout();
}


MagSpikeBox::MagSpikeBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(mod, title, pos, size, "OXYNEURO")
{
	int labelwidth = 60;
	column = 0;
	boxtag = "OXYNEURO";

	InitMenu();

	//SetModFlag(ID_dendstore, "dvstoreflag", "DV Store", 0); 
	//SetModFlag(ID_recep, "recepflag", "Dynamic Receptors", 0); 

	SetModFlag(ID_ipInfusion, "ipInfusionflag", "ip Infus", 0); 
	SetModFlag(ID_ivInfusion, "ivInfusionflag", "iv Infus", 0); 
	SetModFlag(ID_epspsynch, "epspsynchflag", "NMDA Synch", 0); 
	SetModFlag(ID_AHP2mode, "AHP2mode", "Ca Thresh AHP2", 0); 
	SetModFlag(ID_dynostore, "dynostoreflag", "Dyno Store", 0); 

	//OxyPanel();
	VasoPanel();

	/*
	paramset.AddCon("NaDiffhalflife", "NaDif hl", 190, 5, 1); // Half life for the diffusion of NaCl between plasma and EVF
	paramset.AddCon("Osmosishalflife", "Osmo hl", 4.3, 0.1, 2); // Half life for the osmotic exchange between ICF and EVF

	paramset.AddCon("NaClInfused", "NaCl(ml*M)", 3.097, 0.5, 3); // Amount of NaCl infused in ml*M  3.097 = 180.96mg over 30 min
	paramset.AddCon("NaClTimeOfInf", "Inf time", 5, 60, 0); // When does the infusion start
	paramset.AddCon("NaClDurationOfInf", "Inf dura", 1800, 60, 0); // How much does the infusion last
	paramset.AddCon("BasalNaConc", "Basal[Na+]", 155, 1, 1); // Basal concentration of Na+ in every compartment.

	paramset.AddCon("OsmoTimeIv", "Osmo t", 10000, 500, 0); // time of the injection
	paramset.AddCon("OsmoShiftDuration", "OsmoDur", 900, 60, 0); // 15 minutes, 900 s by default
	paramset.AddCon("OsmoNew", "OsmoNew", 305, 1, 0); // New Osmolality induced by hypertonic i.v. or water loading
	paramset.AddCon("OsmoMean", "OsPSP Av", 1, 1, 1); // Mean change of PSP
	paramset.AddCon("OsmoSD", "OsPSP SD", 1, 1, 4); // SD change of PSP
	paramset.AddCon("Osmohalflife", "OsmoHL", 10000, 1000, 0); // Osmotic pressure half life (~ 3 hours). 
	paramset.AddCon("Osmobicucu", "Bicucu", 0, 1, 0); // Bicuculine protocol.
	*/

	paramset.AddCon("synvar", "Syn Var", 1, 0.001, 4);

	ParamLayout(2);

	buttonbox = new wxBoxSizer(wxHORIZONTAL);

	SetPanel(ID_Sec, mod->secbox);
	AddButton(ID_Sec, "Sec", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	SetPanel(ID_Dend, mod->dendbox);
	AddButton(ID_Dend, "Dend", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	SetPanel(ID_Synth, mod->synthbox);
	AddButton(ID_Synth, "Synth", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	SetPanel(ID_Generate, mod->genbox);
	AddButton(ID_Generate, "Gen", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	mainbox->AddSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer();
	mainbox->Add(buttonbox, 0, wxALIGN_CENTRE_HORIZONTAL | wxALIGN_CENTRE_VERTICAL | wxALL, 0);
	mainbox->AddSpacer(5);

	panel->Layout();
}


void MagSpikeBox::OxyPanel()
{
	paramset.AddCon("hstep", "hstep", 1, 1, 2);
	paramset.AddCon("Vrest", "Vrest", -56, 0.1, 4);  // -62
	paramset.AddCon("Vthresh", "Vthresh", -50, 0.1, 4);
	paramset.AddCon("pspmag", "PSP mag", 2, 0.1, 4);  // 1
	paramset.AddCon("psprate", "PSP rate", 190, 1, 2);
	paramset.AddCon("iratio", "Iratio", 0.5, 0.1, 5);
	paramset.AddCon("halflifeMem", "Mem HL", 3.5, 0.1, 4);  // |Like in the Endocrinology paper
	paramset.AddCon("kHAP", "HAP k", 30, 1, 3);
	paramset.AddCon("halflifeHAP", "HAP HL", 7.5, 1, 2); // |Like in the Endocrinology paper
	paramset.AddCon("kDAP", "DAP k", 0, 1, 3);
	paramset.AddCon("halflifeDAP", "DAP HL", 150, 1, 2);
	paramset.AddCon("kAHP", "AHP k", 1, 0.1, 3);
	paramset.AddCon("halflifeAHP", "AHP HL", 350, 1, 2); // Like in the Endocrinology paper

	paramset.AddCon("pspmag2", "PSP2 mag", 0, 0.1, 4); 
	paramset.AddCon("psprate2", "PSP2 rate", 0, 1, 2);
	paramset.AddCon("halflifePSP2", "PSP2 HL", 5, 0.1, 4);

	(*modflags)["modmode"] = 0;
}


void MagSpikeBox::VasoPanel()
{
	paramset.AddCon("hstep", "hstep", 1, 1, 2);
	paramset.AddCon("Vrest", "Vrest", -56, 0.1, 4);  
	paramset.AddCon("Vthresh", "Vthresh", -50, 0.1, 4);
	paramset.AddCon("pspmag", "PSP mag", 2, 0.1, 4); 
	paramset.AddCon("psprate", "PSP rate", 190, 1, 2);
	paramset.AddCon("iratio", "Iratio", 0.5, 0.1, 5);
	paramset.AddCon("halflifeMem", "Mem HL", 3.5, 0.1, 4); 
	paramset.AddCon("kHAP", "HAP k", 30, 1, 3);
	paramset.AddCon("halflifeHAP", "HAP HL", 7.5, 1, 2); 
	paramset.AddCon("kDAP", "DAP k", 0, 1, 3);
	paramset.AddCon("halflifeDAP", "DAP HL", 150, 1, 2);
	paramset.AddCon("kAHP", "AHP k", 1, 0.1, 3);
	paramset.AddCon("halflifeAHP", "AHP HL", 350, 1, 2); 
	paramset.AddCon("kAHP2", "AHP2 k", 1, 0.1, 5);
	paramset.AddCon("halflifeAHP2", "AHP2 HL", 350, 1, 2); 
	paramset.AddCon("aAHP2", "AHP2 a", 0, 1, 5);
	paramset.AddCon("kDyno", "Dyno k", 1.7, 0.1, 4);
	paramset.AddCon("halflifeDyno", "Dyno HL", 10000, 10, 0);
	paramset.AddCon("kCa", "Ca k", 10, 0.1, 2);
	paramset.AddCon("halflifeCa", "Ca HL", 2500, 1, 0);
	paramset.AddCon("gKL", "G Kleak", 16, 0.01, 2);
	paramset.AddCon("gOsmo", "g Osmo", 0, 0.001, 4);
	paramset.AddCon("Ca_rest", "Ca rest", 113, 1, 2);
	paramset.AddCon("ka", "ka", 36, 1, 2);

	paramset.AddCon("pspmag2", "PSP2 mag", 0, 0.1, 4); 
	paramset.AddCon("psprate2", "PSP2 rate", 0, 1, 2);
	paramset.AddCon("halflifePSP2", "PSP2 HL", 5, 0.1, 4);

	(*modflags)["modmode"] = 1;
}


// Copy ParamStore object to panel parameter controls   - moved to ParamBox for general use 27/7/20
// depends on matching tag strings
/*void MagSpikeBox::NeuroParams(ParamStore *params)
{
	int i, conref;
	wxString tag;
	double pval;

	for(i=0; i<paramset.numparams; i++) {
		tag = paramset.con[i]->name;
		pval = (*params)[tag];
		paramset.con[i]->SetValue(pval);
	}
}*/


MagDendBox::MagDendBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(mod, title, pos, size, "OXYDEND")
{
	int labelwidth = 60;
	column = 0;
	boxtag = "OXYDEND";

	InitMenu();

	paramset.AddCon("kstoreDyno", "k Dyn Store", 0.1, 0.01, 5, labelwidth);
	paramset.AddCon("halflifestoreDyno", "HL Dyn Store", 10000, 100, 0, labelwidth);
	paramset.AddCon("spikeDyno", "spike Dyno", 0.01, 0.001, 4, labelwidth);
	paramset.AddCon("tauDynoup", "Dyn tau up", 1, 1, 0, labelwidth);

	paramset.AddCon("kdendCa", "k Dend Ca", 0.1, 0.01, 5, labelwidth);
	paramset.AddCon("halflifedendCa", "HL Dend Ca", 10000, 100, 0, labelwidth);

	ParamLayout(2);

	mainbox->AddStretchSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer(10);
	panel->Layout();

	Connect(wxEVT_CLOSE_WINDOW, wxCloseEventHandler(MagDendBox::OnClose));
}


MagSecBox::MagSecBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(mod, title, pos, size, "OXYSEC")
{
	int labelwidth = 60;
	column = 0;
	boxtag = "OXYSEC";

	InitMenu();

	SetModFlag(ID_diffusion, "diff_flag", "Diffusion", 1); 
	SetModFlag(ID_secfix, "secfix", "Fixed secretion rate", 0);
	//SetModFlag(ID_recep, "recepflag", "Dynamic Receptors", 0); 

	// Secretion - Jorge labels
	/*
	paramset.AddCon("kB", "Broadng", 0.021, 0.001, 3); 
	paramset.AddCon("halflifeB", "Br HL", 2000, 50, 0); 
	paramset.AddCon("Bbase", "Bbasal", 0.5, 0.05, 2);
	paramset.AddCon("kC", "SlowCa K", 0.0003, 0.00001, 5); 
	paramset.AddCon("halflifeC", "C HL(s)", 20000, 1000, 0); 
	paramset.AddCon("kE", "FastCa K", 1.5, 0.02, 2);
	paramset.AddCon("halflifeE", "E HL (ms)", 100, 5, 1); 
	paramset.AddCon("Cth", "CTh", 0.14, 0.01, 3);  
	paramset.AddCon("Cgradient", "C Grad", 5, 0.1, 2);
	paramset.AddCon("Eth", "ETh", 12, 0.05, 2); 
	paramset.AddCon("Egradient", "E Grad", 5, 0.1, 2); 
	paramset.AddCon("beta", "beta", 120, 1, 0.1);
	paramset.AddCon("Rmax", "Res Max", 2000000, 100000, 0);  
	paramset.AddCon("Pmax", "Pool Max", 5000, 500, 0);  
	paramset.AddCon("alpha", "alpha", 0.003, 0.0001, 5); 
	paramset.AddCon("DiffHL", "Ox Diff HL", 61, 5, 0); // 100sec, half life to pass between plasma and ECF. Just a guess.
	paramset.AddCon("ClearHL", "Ox Clear HL", 68, 5, 0); // 58sec half life to be destroyed through the kidneys.
	paramset.AddCon("VolPlasma", "Plasma (ml)", 8.5, 0.5, 1); // Total amount of plasma in a rat. 8.5ml for a 250g rat. 
	paramset.AddCon("VolEVF", "EVFluid (ml)", 9.75, 0.5, 2); // Total amount of Extra Cellular Fluid (without plasma) in a rat. From Fabian et. al (1969) VD = 7.3ml/100g
	paramset.AddCon("secExp", "Sec Exp", 2, 0.1, 2);  // Exponent of the fast [Ca2+], e, when calculating the final secretion.*/

	// Secretion
	paramset.AddCon("kB", "kB", 0.021, 0.001, 3); 
	paramset.AddCon("halflifeB", "halflifeB", 2000, 50, 0); 
	paramset.AddCon("Bbase", "Bbase", 0.5, 0.05, 2);
	paramset.AddCon("kC", "SlowCa K", 0.0003, 0.00001, 5); 
	paramset.AddCon("halflifeC", "C HL", 20000, 1000, 0); 
	paramset.AddCon("kE", "FastCa K", 1.5, 0.02, 2);
	paramset.AddCon("halflifeE", "E HL", 100, 5, 1); 
	paramset.AddCon("Cth", "Cth", 0.14, 0.01, 3);  
	paramset.AddCon("Cgradient", "C Grad", 5, 0.1, 2);
	paramset.AddCon("Eth", "Eth", 12, 0.05, 2); 
	paramset.AddCon("Egradient", "E Grad", 5, 0.1, 2); 
	paramset.AddCon("beta", "beta", 120, 1, 0.1);
	paramset.AddCon("Rmax", "Res Max", 2000000, 100000, 0); 
	paramset.AddCon("Rinit", "Res Init", 2000000, 100000, 0);
	paramset.AddCon("Pmax", "Pool Max", 5000, 500, 0);  
	paramset.AddCon("alpha", "alpha", 0.003, 0.0001, 6); 
	paramset.AddCon("plasma_hstep", "hstep Plas", 1, 1, 0); 
	paramset.AddCon("DiffHL", "Oxy Diff HL", 61, 5, 0); // 100sec, half life to pass between plasma and ECF. Just a guess.
	paramset.AddCon("ClearHL", "Oxy Clear HL", 68, 5, 0); // 58sec half life to be destroyed through the kidneys.
	paramset.AddCon("VolPlasma", "Plasma (ml)", 8.5, 0.5, 1); // Total amount of plasma in a rat. 8.5ml for a 250g rat. 
	paramset.AddCon("VolEVF", "EVFluid (ml)", 9.75, 0.5, 2); // Total amount of Extra Cellular Fluid (without plasma) in a rat. From Fabian et. al (1969) VD = 7.3ml/100g
	paramset.AddCon("secExp", "Sec Exp", 2, 0.1, 2);  // Exponent of the fast [Ca2+], e, when calculating the final secretion.
	paramset.AddCon("secXfix", "secXfix", 0, 0.001, 5);

	ParamLayout(2);

	//wxBoxSizer *paramfilebox = StoreBox("test1");

	mainbox->AddSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer();
	//mainbox->AddSpacer(5);
	//mainbox->Add(paramfilebox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	//mainbox->AddSpacer(2);

	panel->Layout();
}


// Box for showing the analysis of a particular neuron of the Network
MagNeuroDataBox::MagNeuroDataBox(MagNetModel *model, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(model, title, pos, size, "OXYNEURODATA")
{
	int datwidth, labelwidth;
	column = 0;
	mod = model;

	neurodex = 0;  
	neurocount = 0;

	// Neuron selection

	datwidth = 50;
	labelwidth = 70;
	label = NumPanel(datwidth, wxALIGN_CENTRE);
	spikes = NumPanel(datwidth, wxALIGN_RIGHT);
	mean = NumPanel(datwidth, wxALIGN_RIGHT);
	freq = NumPanel(datwidth, wxALIGN_RIGHT);
	sd = NumPanel(datwidth, wxALIGN_RIGHT);

	wxGridSizer *datagrid = new wxGridSizer(2, 5, 5);
	datagrid->Add(new wxStaticText(activepanel, -1, "Name"), 0, wxALIGN_CENTRE);
	datagrid->Add(label);
	datagrid->Add(new wxStaticText(activepanel, -1, "Spikes"), 0, wxALIGN_CENTRE);
	datagrid->Add(spikes);
	datagrid->Add(new wxStaticText(activepanel, -1, "Freq"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(freq);
	datagrid->Add(new wxStaticText(activepanel, -1, "Mean"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(mean);
	datagrid->Add(new wxStaticText(activepanel, -1, "Std Dev"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(sd);

	datneuron = new wxTextCtrl(activepanel, ID_Neuron, "---", wxDefaultPosition, wxSize(50, -1), wxALIGN_LEFT|wxBORDER_SUNKEN|wxST_NO_AUTORESIZE|wxTE_PROCESS_ENTER);
	datspin = new wxSpinButton(activepanel, wxID_ANY, wxDefaultPosition, wxSize(40, 17), wxSP_HORIZONTAL|wxSP_ARROW_KEYS);
	datspin->SetRange(-1000000, 1000000);
	wxBoxSizer *datbox = new wxBoxSizer(wxHORIZONTAL);
	datbox->Add(datspin, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	datbox->AddSpacer(5);

	wxBoxSizer *neurobox = new wxBoxSizer(wxHORIZONTAL);
	neurobox->Add(new wxStaticText(activepanel, wxID_ANY, "Neuron"), 1, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	neurobox->Add(datneuron, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);

	wxStaticBoxSizer *databox = new wxStaticBoxSizer(wxVERTICAL, activepanel, "");
	databox->AddSpacer(2);
	databox->Add(neurobox, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL| wxALL, 5);
	databox->AddSpacer(5);
	databox->Add(datbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL| wxALL, 0);
	databox->AddSpacer(5);
	databox->Add(datagrid, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);

	mainbox->AddSpacer(5);
	mainbox->Add(databox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	mainbox->AddSpacer(2);

	panel->Layout();

	Connect(wxEVT_COMMAND_TEXT_ENTER, wxCommandEventHandler(MagNeuroDataBox::OnEnter));
	Connect(wxEVT_SPIN_UP, wxSpinEventHandler(MagNeuroDataBox::OnNext));
	Connect(wxEVT_SPIN_DOWN, wxSpinEventHandler(MagNeuroDataBox::OnPrev));
	Connect(wxEVT_SPIN, wxSpinEventHandler(MagNeuroDataBox::OnSpin));
}


// Jorge comment - returns calculations for each neuron of the network when we want to see their graphs
void MagNeuroDataBox::NeuroData()
{

	// To show results, numerical and graphically, we need to:
	//		- send the spiketimes and number of spikes of the neuron we want the spyke statitistic to neurocalc
	//		- send the secretion and plasma data to single vectors to the graph class (it does not work with arrays of vectors)
	mod->currmodneuron->neurocalc(&(mod->modneurons[neurodex]));
	mod->currmodneuron->id = neurodex;
	//mod->neurons->index = neurodex;    // what is this doing?   25/11/20

	if(mod->burstbox) mod->burstbox->ModDataScan();

	// Assigning to the graphbase vector with name "OxSecretion" the data from the position neurodex of the variable OxSecretion of the class array neurons
	(*mod->graphbase)["Secretion"]->gdatadv = &(mod->modneurons[neurodex].Secretion);  // gdatadv: graph data double vector. 
	//(*mod->graphbase)["OxyPlasma"]->gdatadv = &(mod->neurons[neurodex].OxyPlasma); 

	mod->magpop->storeLong = mod->modneurons[neurodex].storeLong;
	mod->magpop->synthstoreLong = mod->modneurons[neurodex].synthstoreLong;
	mod->magpop->synthrateLong = mod->modneurons[neurodex].synthrateLong;
	mod->magpop->secLong = mod->modneurons[neurodex].secLong;
	mod->magpop->secHour = mod->modneurons[neurodex].secHour;
	mod->magpop->transLong = mod->modneurons[neurodex].transLong;

	PanelData(&(mod->modneurons[neurodex]));	

	mod->spikebox->CopyParams(mod->modneurons[neurodex].spikeparams);
	mod->secbox->CopyParams(mod->modneurons[neurodex].secparams);
	mod->synthbox->CopyParams(mod->modneurons[neurodex].synthparams);
}


void MagNeuroDataBox::PanelData(NeuroDat *data)
{
	wxString snum;

	if(data->netflag) snum = "sum";
	else snum = numstring(neurodex, 0);
	datneuron->SetLabel(snum);

	label->SetLabel(data->name);
	spikes->SetLabel(snum.Format("%d", data->spikecount));
	freq->SetLabel(snum.Format("%.2f", data->freq));
	mean->SetLabel(snum.Format("%.1f", data->meanisi));
	sd->SetLabel(snum.Format("%.2f", data->isivar));

}


void MagNeuroDataBox::OnPrev(wxSpinEvent& WXUNUSED(event))
{
	// wxSpinButton Diagnostic
	//diagbox->Write(text.Format("spin down neurodex %d  val %d  min %d  max %d\n", neurodex, datspin->GetValue(), datspin->GetMin(), datspin->GetMax()));

	if(!neurocount) return;
	if(neurodex > 0) neurodex--;  // going to the previous neuron
	else neurodex = neurocount-1;

	NeuroData();
    mainwin->scalebox->GraphUpdate();
}


void MagNeuroDataBox::OnNext(wxSpinEvent& WXUNUSED(event))
{
	// wxSpinButton Diagnostic
	//diagbox->Write(text.Format("spin up neurodex %d  val %d  min %d  max %d\n", neurodex, datspin->GetValue(), datspin->GetMin(), datspin->GetMax()));

	if(!neurocount) return;
	if(neurodex < neurocount-1) neurodex++;
	else neurodex = 0; // index 0 points to the first neuron

	NeuroData();
    mainwin->scalebox->GraphUpdate();
}


void MagNeuroDataBox::OnEnter(wxCommandEvent& event)
{
	int id = event.GetId();
	long data;

	// Enter pressed for neuron selection (to see what is happenning with a particular neuron, selected by its index)
	if(id == ID_Neuron) {
		datneuron->GetValue().ToLong(&data);
		if(data >= 0 && data < neurocount) {
			neurodex = data;
			NeuroData();
		}
		return;
	}
	else NeuroData();
    mainwin->scalebox->GraphUpdate();
}


// MagNetGridBox - derived version of GridBox class

MagNetGridBox::MagNetGridBox(MagNetModel *model, const wxString& title, const wxPoint& pos, const wxSize& size, int rows, int cols)
	: GridBox(model, title, pos, size, rows, cols, true, true)
{
	//GridDefault();
	//mod = model;
}


void MagNetGridBox::OnPlot(wxCommandEvent& event)
{
	mod->RangePlot(currgrid);
	WriteVDU("Plot\n");
	//diagbox->Write("param scan\n");
}



MagNetProtoBox::MagNetProtoBox(MagNetModel *model, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(model, title, pos, size, "PROTO", 1, 1)
{
	int pnum, inpnum, rampnum;
	int pulsenum0, pulsenum1, rampnum0, rampnum1;
	int numwidth;
	boxtag = "PROTO";
	mod = model;
	wxString tag;

	long notestyle = wxAUI_NB_TOP | wxAUI_NB_TAB_SPLIT | wxAUI_NB_TAB_MOVE | wxAUI_NB_SCROLL_BUTTONS;
	wxAuiNotebook *tabpanel = new wxAuiNotebook(this, wxID_ANY, wxDefaultPosition, wxDefaultSize, notestyle);

	pnum = 0;
	//status = mod->netbox->status;


	// Range Panel

	ToolPanel *rangepanel = new ToolPanel(this, tabpanel);
	rangepanel->SetFont(boxfont);
	wxBoxSizer *rangesizer = new wxBoxSizer(wxVERTICAL);
	rangepanel->SetSizer(rangesizer);

	activepanel = rangepanel;
	paramset.panel = activepanel;

	labelwidth = 50;
	numwidth = 45;

	paramset.AddNum("rangestart", "Start", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rangestop", "Stop", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rangestep", "Step", 300, 0, labelwidth, numwidth); 
	paramset.AddNum("rangedata", "Data", 300, 0, labelwidth, numwidth); 

	wxStaticBoxSizer *rangebox0 = new wxStaticBoxSizer(wxVERTICAL, rangepanel, "Range Input");
	for(pnum=pnum; pnum<paramset.numparams; pnum++) {
		rangebox0->Add(paramset.con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
	}

	wxBoxSizer *rangerunbox = new wxBoxSizer(wxHORIZONTAL);
	rangerunbox->Add(TextLabel("Input"), 1, wxALIGN_CENTRE);
	currentrange = NumPanel(40, wxALIGN_CENTRE, "---"); 
	rangerunbox->AddSpacer(10);
	rangerunbox->Add(currentrange, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	rangebox0->AddSpacer(10);
	rangebox0->Add(rangerunbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);

	rangebox0->AddSpacer(10);
	AddButton(ID_Range, "Run", 50, rangebox0);

	wxBoxSizer *rangebox = new wxBoxSizer(wxHORIZONTAL);
	rangebox->Add(rangebox0, 0, wxALL, 5);
	rangesizer->AddSpacer(10);
	rangesizer->Add(rangebox, 1, wxALL, 0);
	rangepanel->Layout();


	// Pulse Panel 

	ToolPanel *pulsepanel = new ToolPanel(this, tabpanel);
	pulsepanel->SetFont(boxfont);
	wxBoxSizer *pulsesizer = new wxBoxSizer(wxVERTICAL);
	pulsepanel->SetSizer(pulsesizer);

	activepanel = pulsepanel;
	paramset.panel = activepanel;

	labelwidth = 50;
	numwidth = 45;

	paramset.AddNum("pulsebase0", "Base", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("pulsestart0", "Start", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("pulsestop0", "Stop", 300, 0, labelwidth, numwidth); 
	paramset.AddNum("pulseinit0", "Initial", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("pulsehl0", "Halflife", 0.1, 2, labelwidth, numwidth); 

	pulsenum0 = paramset.numparams;

	paramset.GetCon("pulseinit0")->SetMinMax(-10000, 10000);

	wxStaticBoxSizer *pulsebox0 = new wxStaticBoxSizer(wxVERTICAL, pulsepanel, "Pulse Input L1");
	for(pnum=pnum; pnum<paramset.numparams; pnum++) {
		pulsebox0->Add(paramset.con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
	}
	wxBoxSizer *pulserunbox = new wxBoxSizer(wxHORIZONTAL);

	pulserunbox->Add(TextLabel("Input"), 1, wxALIGN_CENTRE);
	currentpulse = NumPanel(40, wxALIGN_CENTRE, "---"); 
	pulserunbox->AddSpacer(10);
	pulserunbox->Add(currentpulse, 1, wxALIGN_CENTRE_VERTICAL);
	pulsebox0->AddSpacer(10);
	pulsebox0->Add(pulserunbox, 1, wxRIGHT|wxLEFT, 5);
	pulsebox0->AddSpacer(10);
	AddButton(ID_Pulse, "Run", 50, pulsebox0);

	paramset.AddNum("pulsebase1", "Base", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("pulsestart1", "Start", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("pulsestop1", "Stop", 300, 0, labelwidth, numwidth); 
	paramset.AddNum("pulseinit1", "Initial", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("pulsehl1", "Halflife", 0.1, 2, labelwidth, numwidth); 

	paramset.GetCon("pulseinit1")->SetMinMax(-10000, 10000);

	wxStaticBoxSizer *pulsebox1 = new wxStaticBoxSizer(wxVERTICAL, pulsepanel, "Pulse Input L2");
	for(pnum=pnum; pnum<paramset.numparams; pnum++) {
		pulsebox1->Add(paramset.con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxRIGHT|wxLEFT, 5);
	}

	pulsebox1->AddSpacer(10);
	AddButton(ID_Pulse, "Run", 50, pulsebox1);


	wxBoxSizer *pulsebox = new wxBoxSizer(wxHORIZONTAL);
	pulsebox->Add(pulsebox0, 0, wxALL, 5);
	pulsebox->AddStretchSpacer();
	pulsebox->Add(pulsebox1, 0, wxALL, 5);
	pulsesizer->AddSpacer(10);
	pulsesizer->Add(pulsebox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	pulsepanel->Layout();


	// Ramp Panel
	ToolPanel *ramppanel = new ToolPanel(this, tabpanel);
	ramppanel->SetFont(boxfont);
	wxBoxSizer *rampsizer = new wxBoxSizer(wxVERTICAL);
	ramppanel->SetSizer(rampsizer);

	activepanel = ramppanel;
	paramset.panel = activepanel;

	labelwidth = 50;
	numwidth = 45;

	paramset.AddNum("rampbase0", "Base", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rampstart0", "Start", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rampstop0", "Stop", 300, 0, labelwidth, numwidth); 
	paramset.AddNum("rampinit0", "Init", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rampafter0", "After", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rampstep0", "1s Step", 0.1, 4, labelwidth, numwidth); 
	paramset.SetMinMax("rampstep0", -1000, 1000);
	paramset.SetMinMax("rampinit0", -1, 100000);
	paramset.SetMinMax("rampafter0", -1, 100000);
	rampnum0 = paramset.numparams;

	wxStaticBoxSizer *rampbox0 = new wxStaticBoxSizer(wxVERTICAL, ramppanel, "Input Ramp L1");
	for(pnum=pnum; pnum<paramset.numparams; pnum++) {
		rampbox0->Add(paramset.con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxRIGHT|wxLEFT, 5);
	}
	wxBoxSizer *inputrunbox = new wxBoxSizer(wxHORIZONTAL);

	inputrunbox->Add(TextLabel("Input"), 1, wxALIGN_CENTRE);
	currentinput = NumPanel(40, wxALIGN_CENTRE, "---"); 
	inputrunbox->AddSpacer(10);
	inputrunbox->Add(currentinput, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	rampbox0->AddSpacer(10);
	rampbox0->Add(inputrunbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
	rampbox0->AddSpacer(10);
	AddButton(ID_Ramp, "Run", 50, rampbox0);

	paramset.AddNum("rampbase1", "Base", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rampstart1", "Start", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rampstop1", "Stop", 300, 0, labelwidth, numwidth); 
	paramset.AddNum("rampmax1", "Max", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("rampgrad1", "Grad", 0.005, 4, labelwidth, numwidth); 

	//paramset.SetMinMax("rampstep1", -1000, 1000);
	rampnum1 = paramset.numparams;

	wxStaticBoxSizer *rampbox1 = new wxStaticBoxSizer(wxVERTICAL, ramppanel, "Input Ramp Curve");
	for(pnum=pnum; pnum<paramset.numparams; pnum++) {
		rampbox1->Add(paramset.con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxRIGHT|wxLEFT, 5);
	}
	rampbox1->AddSpacer(10);
	AddButton(ID_RampCurve, "Run", 50, rampbox1);

	wxBoxSizer *rampbox = new wxBoxSizer(wxHORIZONTAL);
	rampbox->Add(rampbox0, 0, wxALL, 5);
	rampbox->AddStretchSpacer();
	rampbox->Add(rampbox1, 0, wxALL, 5);
	rampsizer->AddSpacer(10);
	rampsizer->Add(rampbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALL, 0);
	ramppanel->Layout();


	// Gavage Panel
	ToolPanel *gavpanel = new ToolPanel(this, tabpanel);
	gavpanel->SetFont(boxfont);
	wxBoxSizer *gavsizer = new wxBoxSizer(wxVERTICAL);
	gavpanel->SetSizer(gavsizer);

	activepanel = gavpanel;
	paramset.panel = activepanel;

	labelwidth = 50;
	numwidth = 45;

	paramset.AddNum("gavbase0", "Base", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("gavstart0", "Start", 200, 0, labelwidth, numwidth); 
	paramset.AddNum("gavstop0", "Stop", 300, 0, labelwidth, numwidth); 
	paramset.AddNum("gavstep0", "1s Step", 0.1, 2, labelwidth, numwidth); 
	//paramset.SetMinMax("gavstep0", -1000, 1000);

	wxStaticBoxSizer *gavbox0 = new wxStaticBoxSizer(wxVERTICAL, gavpanel, "Gavage");
	for(pnum=pnum; pnum<paramset.numparams; pnum++) {
		gavbox0->Add(paramset.con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxRIGHT|wxLEFT, 5);
	}
	wxBoxSizer *gavrunbox = new wxBoxSizer(wxHORIZONTAL);

	//gavrunbox->Add(TextLabel("Input"), 1, wxALIGN_CENTRE);
	//currgavinput = NumPanel(40, wxALIGN_CENTRE, "---"); 
	//gavrunbox->AddSpacer(10);
	//gavrunbox->Add(currentinput, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	gavbox0->AddSpacer(10);
	gavbox0->Add(gavrunbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
	gavbox0->AddSpacer(10);
	AddButton(ID_Gavage, "Run", 50, gavbox0);

	wxBoxSizer *gavbox = new wxBoxSizer(wxHORIZONTAL);
	gavbox->Add(gavbox0, 0, wxALL, 5);
	gavbox->AddStretchSpacer();
	//gavbox->Add(rampbox1, 0, wxALL, 5);
	gavsizer->AddSpacer(10);
	gavsizer->Add(gavbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALL, 0);
	gavpanel->Layout();


	//////////////////////////////////////////////////
	// Main Structure

	tabpanel->Freeze();
	tabpanel->AddPage(ramppanel, "Ramp" , true);
	tabpanel->AddPage(pulsepanel, "Pulse" , false);
	tabpanel->AddPage(rangepanel, "Range" , false);
	tabpanel->AddPage(gavpanel, "Gavage" , false);
	tabpanel->Thaw();

	ToolPanel *storepanel = new ToolPanel(this, wxDefaultPosition, wxDefaultSize);
	wxBoxSizer *storesizer = new wxBoxSizer(wxVERTICAL);
	storepanel->SetSizer(storesizer);

	activepanel = storepanel;
	wxBoxSizer *paramfilebox = StoreBox("ns1", storepanel);

	storesizer->Add(paramfilebox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	storesizer->Layout();

	winman->AddPane(tabpanel, wxAuiPaneInfo().Name("tabpane").CentrePane().PaneBorder(false));
	winman->AddPane(storepanel, wxAuiPaneInfo().Name("storepane").Bottom().CaptionVisible(false).BestSize(-1, 30).PaneBorder(false));
	winman->Update();

	Connect(ID_Ramp, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetProtoBox::OnRun));
	Connect(ID_Pulse, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetProtoBox::OnRun));
	Connect(ID_Range, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetProtoBox::OnRun));
	Connect(ID_Gavage, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetProtoBox::OnRun));
	Connect(ID_RampCurve, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagNetProtoBox::OnRun));
}


void MagNetProtoBox::OnRun(wxCommandEvent& event)
{
	if(event.GetId() == ID_Ramp) (*mod->modeflags)["prototype"] = ramp;
	if(event.GetId() == ID_Pulse) (*mod->modeflags)["prototype"] = pulse;
	if(event.GetId() == ID_Range) (*mod->modeflags)["prototype"] = range;
	if(event.GetId() == ID_Gavage) (*mod->modeflags)["prototype"] = gavage;
	if(event.GetId() == ID_RampCurve) (*mod->modeflags)["prototype"] = rampcurve;

	//mod->netbox->SetNeuroCount();
	mod->netbox->countmark = 0;
	mod->RunModel();
} 

// Control box for testing signal processing
MagSignalBox::MagSignalBox(MagNetModel *vmnmodel, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(vmnmodel, title, pos, size, "Signal")
{
	column = 0;
	boxtag = "Signal";
	mod = vmnmodel;

	InitMenu();

	//SetModFlag(ID_revpots, "revpots", "Reversal Potentials", 0); 
	SetModFlag(ID_noise, "noiseflag", "Noise Signal", 0); 


	// Parameter controls
	//
	// AddCon(tag string, display string, initial value, click increment, decimal places)
	// ----------------------------------------------------------------------------------


	paramset.AddCon("noimean", "Noise Mean", 300, 1, 2);
	paramset.AddCon("noitau", "Noise Tau", 1000, 1, 2);
	paramset.AddCon("noiamp", "Noise Amp", 1, 0.1, 2);
	paramset.AddCon("sigIratio", "Iratio", 0, 0.1, 2);

	/*if(!mod->basicmode) {
	paramset.AddCon("synwaveamp", "Wave Amp", 0, 1, 2);
	paramset.AddCon("synwavecycle", "Wave Cycle", 1000, 1, 2);
	paramset.AddCon("synwaveshift", "Wave Shift", 0, 0.1, 2);

	paramset.AddCon("kB", "kB", 0, 1, 4);
	paramset.AddCon("halflifeB", "halflifeB", 100, 1, 2);
	paramset.AddCon("timerange", "timerange", 100, 1, 0);
	}*/

	ParamLayout(2);

	// ----------------------------------------------------------------------------------

	defbutt = 0;
	wxBoxSizer *runbox = RunBox();

	//wxSizer *paramfilebox = StoreBox("test1");

	mainbox->AddSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer();
	mainbox->Add(runbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	//mainbox->AddSpacer(5);
	//mainbox->Add(paramfilebox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	//mainbox->AddSpacer(2);

	panel->Layout();
}


void MagSignalBox::OnRun(wxCommandEvent& event)
{
	//mod->SigSim(mod->netdat1);
	//mod->SigSim(mod->netdat2);

	//ParamStore *sigparams = GetParams();
	//int timerange = (*sigparams)["timerange"];
	//if((*modflags)["layerflag"] == 0) mod->currvmn->MeanSpikeForm(mod->neurodata->V, timerange, (*mod->netbox->modflags)["formfilter"]);
	//if((*modflags)["layerflag"] == 1) mod->currvmn->MeanSpikeForm(mod->neurodata->V2, timerange, (*mod->netbox->modflags)["formfilter"]);

	mod->RunModel();
	mainwin->scalebox->GraphUpdate();
}



MagGenBox::MagGenBox(MagNetModel *model, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(model, title, pos, size)
{
	labelwidth = 50;
	int sdwidth = 30;
	boxtag = "MAGNOGEN";
	mod = model;

	//SetMenuBar(menuBar);

	//gentags[0] = "halflifeHAP";
	//gentags[1] = "kAHP";
	//gentags[2] = "kAHP2";
	//gentags[3] = "kDAP";
	//gentags[4] = "halflifeDyno";
	//gentags[5] = "ratioDyno";
	//gentags[6] = "kCa";
	//gentags[7] = "gKL";
	//gentags[8] = "pspmag";
	//gentags[9] = "Vrest";
	//gentags[10] = "synvar";

	numgen = 2;

	//paramset.AddCon("halflifeHAPbase", "HAP hl", 9, 1, 1, labelwidth); 
	//paramset.AddCon("halflifeHAPsd", "SD", 1, 0.1, 1, sdwidth); 

	gentags[0] = "kAHP";
	paramset.AddCon("kAHPbase", "kAHP", 0.04, 0.005, 5, labelwidth);
	paramset.AddCon("kAHPsd", "SD", 0.01, 0.001, 5, sdwidth);

	//paramset.AddCon("kAHP2base", "kAHP2", 0.00001, 0.005, 6, labelwidth);
	//paramset.AddCon("kAHP2sd", "SD", 0.00001, 0.001, 6, sdwidth);

	gentags[1] = "kDAP";
	paramset.AddCon("kDAPbase", "kDAP", 1, 0.1, 2, labelwidth);
	paramset.AddCon("kDAPsd", "SD", 0.5, 0.1, 2, sdwidth);

	/*
	paramset.AddCon("ratioDynobase", "Dyno rat", 14, 0.1, 2, labelwidth);
	paramset.AddCon("ratioDynosd", "SD", 1.5, 0.1, 2, sdwidth);

	paramset.AddCon("halflifeDynobase", "Dyno hl", 7500, 500, 0, labelwidth);
	paramset.AddCon("halflifeDynosd", "SD", 1000, 100, 0, sdwidth);

	paramset.AddCon("kCabase", "kCa", 11, 0.5, 1, labelwidth);
	paramset.AddCon("kCasd", "SD", 1, 0.5, 1, sdwidth);

	paramset.AddCon("gKLbase", "gKL", 17, 1, 1, labelwidth);
	paramset.AddCon("gKLsd", "SD", 2, 0.5, 1, sdwidth);

	paramset.AddCon("pspmagbase", "pspmag", 4, 0.1, 2, labelwidth);
	paramset.AddCon("pspmagsd", "SD", 0, 0.1, 2, sdwidth);

	paramset.AddCon("Vrestbase", "Vrest", -62, 0.1, 2, labelwidth);
	paramset.AddCon("Vrestsd", "SD", 0, 0.1, 2, sdwidth);

	paramset.AddCon("synvarbase", "synvar", 1, 0.1, 2, labelwidth);
	paramset.AddCon("synvarsd", "SD", 0, 0.05, 2, sdwidth);
	*/


	wxFlexGridSizer *gengrid = new wxFlexGridSizer(2, 5, 0);
	for(i=0; i<numgen; i++) {
		//gengrid->Add(TextLabel(labelset[i]), 0, wxALIGN_CENTRE);
		gengrid->Add(paramset.con[i*2]);
		gengrid->Add(paramset.con[i*2+1]);
	}

	parambox->AddSpacer(10);
	parambox->Add(gengrid, 0);

	buttonbox = new wxBoxSizer(wxHORIZONTAL);
	AddButton(ID_Generate, "Generate", 80, buttonbox);
	AddButton(ID_Zero, "Zero", 80, buttonbox);

	mainbox->AddSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddSpacer(10);
	mainbox->Add(buttonbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);	

	panel->Layout();

	Connect(ID_Zero, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagGenBox::OnZero));
	Connect(ID_Generate, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagGenBox::OnGenerate));
	Connect(wxEVT_CLOSE_WINDOW, wxCloseEventHandler(MagGenBox::OnClose));
	Connect(ID_paramstore, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagGenBox::OnParamStore));
	Connect(ID_paramload, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MagGenBox::OnParamLoad));
}


void MagGenBox::OnZero(wxCommandEvent& WXUNUSED(event))
{
	for(i=0; i<numgen; i++)
		paramset.con[(int)paramset.ref[gentags[i] + "sd"]]->SetValue(0);
}


void MagGenBox::OnGenerate(wxCommandEvent& WXUNUSED(event))
{
	ParamStore *genparams = GetParams();
	//mod->spikebox->NeuroGen(genparams);
}


void MagGenBox::OnClose(wxCloseEvent& event)
{
	this->Show(false);
}