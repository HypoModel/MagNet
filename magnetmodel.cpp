#include "magnetmodel.h"



MagNetDat::MagNetDat()
{
}


MagNetModel::MagNetModel(int type, wxString name, HypoMain *main)
: NeuroMod(type, name, main)
{
	int i;

	path = "MagNet";
	oldhist = false;
    projmode = true;

	storesize = 100000;
	diagbox = mainwin->diagbox;

	mainwin->SetMinSize(wxSize(470, 300));
	mainwin->xstretch = 0;
	xmin = -1000000;

	datsample = 1000;    // recording rate
	popscale = 9000;     // secretion scaling from single neuron to in vivo population

	netready = false;
	modneurons.resize(200);
	modneurons_max = 200;

	currmodneuron = new SpikeDat();
	currmodneuron->BurstInit();
	neuroindex = 0;

	neurodata = new MagNeuroDat();

	netdata = new MagNetDat();
	magpop = new MagPop();
	netdat = new SpikeDat();
	netneuron = new SpikeDat();

	// Diagnostic long spike record
	//neurons[0].diagstore();
	//currmodneuron->maxspikes = 1000000;
	//currmodneuron->times.resize(1000000);
	//netneuron->maxspikes = 1000000;
	//netneuron->times.resize(1000000);

	// Grid Plot Data
	for(i=0; i<10; i++) {
		gridplot[i].setsize(1000);
		gridplotx[i].setsize(1000);
		gridploterr[i].setsize(1000);
	}

	// Protocol Data
	rangecount = 10;
	rangedatamax = 1000;
	for(i=0; i<rangecount; i++) rangedata[i].setsize(rangedatamax);
	rangeref.setsize(rangedatamax);   // common X values for plotting range protocol data
	rangeindex = 0;

	celltypes = 1;    // reserve for possible future multiple cell types as in VMNNet

	// Initialise neuro select data
	magpop->OxySecretion = &modneurons[0].Secretion;
	magpop->OxyPlasma = &modneurons[0].Plasma;
	
	gridbox = new MagNetGridBox(this, "Data Grid", wxPoint(0, 0), wxSize(320, 500), 100, 20);
	neurobox = new NeuroBox(this, "Spike Data", wxPoint(0, 0), wxSize(320, 500));
	secbox = new MagSecBox(this, "Secretion and Diffusion", wxPoint(0, 0), wxSize(320, 500));
	dendbox = new MagDendBox(this, "Dendritic", wxPoint(0, 0), wxSize(320, 500));
	neurodatabox = new MagNeuroDataBox(this, "Model Neuron Data", wxPoint(0, 0), wxSize(320, 500));
	signalbox = new MagSignalBox(this, "Signal Box", wxPoint(0, 300), wxSize(400, 500));
	protobox = new MagNetProtoBox(this, "Protocol", wxPoint(0, 0), wxSize(320, 500));
	synthbox = new MagSynthBox(this, "Synthesis", wxPoint(0, 0), wxSize(320, 500));

	// Panel control boxes, must come last to link panel buttons
	spikebox = new MagSpikeBox(this, "Spiking", wxPoint(0, 0), wxSize(320, 500));
	netbox = new MagNetBox(this, mainwin, "Hypo Net Model", wxPoint(0, 0), wxSize(320, 500));

	// link mod owned boxes
	mainwin->neurobox = neurobox;
	mainwin->gridbox = gridbox;

	// general HypoMod modules
	mainwin->BurstModule(this, currmodneuron); 
	burstbox = mainwin->burstbox;

	mainwin->PlotModule(this); 
	plotbox = mainwin->plotbox;

	//oxynetbox->canclose = false;
	dispbox = netbox;

	// Cell data NeuroDat store vector and view SpikeDat
	int numview = 2;
	viewcell.resize(numview); 
	viewcell[0].BurstInit();
	celldata.resize(10);

	// NeuroBox panel for multi cell data import, analysis and fitting
	neurobox->textgrid = gridbox->textgrid[0];
	neurobox->gridbox = gridbox;
	neurobox->burstbox = mainwin->burstbox;

	// Cell data linking
	neurobox->cellpanel->SetData(&viewcell[0], &celldata);

	// Mod data linking
	neurobox->AddModSpikePanel(currmodneuron, (std::vector<NeuroDat>*)&modneurons);

	// GridBox linking and set up
	gridbox->celldata = &celldata;
	gridbox->neurobox = neurobox;

	gridbox->NeuroButton();
	gridbox->PlotButton();
	gridbox->ParamButton();


#ifdef HYPOSOUND
	mainwin->SoundModule(this);
	soundbox = mainwin->soundbox;
#endif

	//(*toolflags)["spikebox"] = 1;		// Select universal model tools
	//mainwin->ToolLoad(this);				// Load universal model tools

	modtools.AddBox(netbox, true);
	modtools.AddBox(spikebox, true);
	modtools.AddBox(secbox, true);
	modtools.AddBox(dendbox, true);
	modtools.AddBox(neurodatabox, true);
	modtools.AddBox(gridbox, true);
	modtools.AddBox(neurobox, true);
	modtools.AddBox(signalbox, true);
	modtools.AddBox(protobox, true);
	modtools.AddBox(synthbox, true);
    #ifdef HYPOSOUND
    modtools.AddBox(soundbox, true);
    #endif
	modbox = netbox;

	//ModLoad();
	/*
	for(i=0; i<modtools.numtools; i++) {
		modtools.box[i]->ReSize();
		modtools.box[i]->Show(modtools.box[i]->visible);
	}*/
	
	netbox->ParamLoad("default");
	graphload = false;

	gsync = 0;

	//oxynetbox->Show(true);
	GraphData();

    Connect(wxID_ANY, EVT_MODTHREAD_COMPLETED, wxThreadEventHandler(MagNetModel::OnModThreadCompletion));   
}


// ParamScan() reads model parameters from textgrid into parambox, parameters across columns
// tags in row 1, values in row 2

void MagNetModel::ParamScan()
{
	int i, paramcount;
	TextGrid *paramgrid;
	wxString text, paramtag;
	double paramval;
	ParamCon *paramcon;

	// Link panel data from GridBox
	if(gridbox->paramgrid) paramgrid = gridbox->paramgrid;
	else {
		gridbox->WriteVDU("No param data found\n");
		return;
	}

	paramcount = 0;

	// reads first parameter tag
	paramtag = paramgrid->GetCell(1, 0);
	paramtag.Trim();

	while(!paramtag.IsEmpty()) {
		paramcon = spikebox->paramset.GetCon(paramtag);
		// check valid parameter tag
		if(!paramcon) diagbox->Write(text.Format("Param %s tag not found\n", paramtag));	
		else {
			// check and read valid parameter value
			if(!paramgrid->CheckDouble(paramgrid->selectrow, paramcount, &paramval)) diagbox->Write(text.Format("Param %s bad param value\n", paramtag));
			else paramcon->SetValue(paramval);
		}
		// read next parameter tag
		paramcount++;
		paramtag = paramgrid->GetCell(1, paramcount);
		paramtag.Trim();		
	}
	gridbox->WriteVDU(text.Format("Parameters read OK\n"));
}


void MagNetModel::OnModThreadCompletion(wxThreadEvent&)
{
    runmute->Lock();
    runflag = 0;
    runmute->Unlock();
    mainwin->scalebox->GraphUpdate();
    
    diagbox->Write("Model thread OK\n\n");
}


void MagNetModel::OnEndRun(wxCommandEvent&)
{
    if(mainwin->diagnostic) diagbox->Write("Model Run Finished\n");
    //mainwin->scalebox->GraphUpdate();
    //runflag = 0;
}


void MagNetModel::ModClose()
{
	mainwin->burstbox->Store();
	neurobox->Store();
#ifdef HYPOSOUND
	soundbox->Store();
#endif
}


void MagNetModel::SoundOn()
{
	(*graphwin)[0].Highlight(50);
}


int MagNetModel::SoundLink(SoundBox *soundbox)
{
#ifdef HYPOSOUND
	if(currmodneuron->spikecount) {
		soundbox->SetSpikeData(currmodneuron);
		soundbox->SetWaveData(&neurodata->Ca);
	}
	else if(viewcell[0].spikecount) {
		soundbox->SetSpikeData(&viewcell[0]);	
	}
	
	return soundbox->spikecount;
#endif
}


void MagNetModel::RangePlot(TextGrid *textgrid)
{
	int i, count;
	wxString text;

	count = gridbox->ColumnData(0, &rangeref);
	for(i=0; i<rangecount; i++) {
		gridbox->ColumnData(i+1, &rangedata[i]);
		graphbase->GetGraph(text.Format("rangedata%d", i))->xcount = count;
	}

	mainwin->scalebox->GraphUpdate();
}


void MagNetModel::ScaleConsoleBelow(ScaleBox *scalebox, int condex)
{
	int ostype = scalebox->ostype;
	wxBoxSizer *vbox = scalebox->consolebox[condex];
	//ToolButton *spikesbutton, *ratebutton;

	scalebox->buttonheight = 20;

	if(condex == 0) {
		wxBoxSizer *resbox = new wxBoxSizer(wxHORIZONTAL); 
		wxBoxSizer *modebox = new wxBoxSizer(wxHORIZONTAL); 

		if(ostype == Mac) {
			//scalebox->ScaleButton(ID_spikes, "Spike", 40, resbox);
			scalebox->GraphButton("spikeres", 0, ID_spikes, "Spikes", 40, resbox); 
			//scalebox->ScaleButton(ID_rateres, "Rate", 40, resbox);
			scalebox->GraphButton("rateres", 0, ID_rateres, "Rate", 40, resbox)->linkID = ID_spikes; 
			//scalebox->GraphButton("nettog", 0, ID_net, "Net", 43, modebox);
			scalebox->databutton = scalebox->ScaleButton(ID_data, "Disp", 43, modebox);
		}
		else {
			//scalebox->ScaleButton(ID_spikes, "Spikes", 37, resbox); 
			scalebox->GraphButton("spikeres", 0, ID_spikes, "Spikes", 37, resbox); 
			resbox->AddSpacer(2);
			//scalebox->ScaleButton(ID_rateres, "Rate", 37, resbox); 
			scalebox->GraphButton("rateres", 0, ID_rateres, "Rate", 37, resbox, 3)->linkID = ID_spikes; 
			//scalebox->GraphButton("nettog", 0, ID_net, "Net", 37, modebox);
			scalebox->databutton = scalebox->ScaleButton(ID_data, "All", 43, modebox);
		}
		vbox->Add(resbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
		vbox->Add(modebox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	}
}


void MagNetModel::ScaleConsoleAbove(ScaleBox *scalebox, int condex)
{
	int ostype = scalebox->ostype;
	wxBoxSizer *vbox = scalebox->consolebox[condex];
    
    scalebox->buttonheight = 20;

	/*if(condex == 0) {
		wxBoxSizer *resbox = new wxBoxSizer(wxHORIZONTAL); 
		wxBoxSizer *modebox = new wxBoxSizer(wxHORIZONTAL); 

		if(ostype == Mac) {
			scalebox->ScaleButton(ID_spikes, "Spike", 40, resbox);
			scalebox->ScaleButton(ID_rateres, "Rate", 40, resbox);
			scalebox->GraphButton("nettog", 0, ID_net, "Net", 43, modebox);
		}
		else {
			scalebox->ScaleButton(ID_spikes, "Spikes", 37, resbox); 
			resbox->AddSpacer(2);
			scalebox->ScaleButton(ID_rateres, "Rate", 37, resbox); 
			scalebox->GraphButton("nettog", 0, ID_net, "Net", 37, modebox);
		}
		vbox->Add(resbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
		vbox->Add(modebox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	}*/


	if(condex == 2) {
		wxBoxSizer *binbox = new wxBoxSizer(wxHORIZONTAL); 
		wxBoxSizer *filtbox = new wxBoxSizer(wxHORIZONTAL); 
		if(ostype == Mac) {
			scalebox->GraphButton("hazmode1", 0, ID_histhaz1, "Hist / Haz", 70, vbox);
			scalebox->GraphButton("binrestog1", 0, ID_binres1, "Bin", 45, binbox);
			scalebox->GraphButton("normtog", 0, ID_norm, "Norm", 45, binbox);
			scalebox->burstbutton = scalebox->GraphButton("burstmode", 0, ID_burst, "Burst", 45, filtbox);
			scalebox->selectbutton = scalebox->GraphButton("selectmode", 0, ID_select, "Select", 45, filtbox);
		}
		else {
			scalebox->GraphButton("hazmode1", 0, ID_histhaz1, "Hist / Haz", 54, vbox);
			scalebox->GraphButton("binrestog1", 0, ID_binres1, "Bin Res", 43, binbox);
			scalebox->GraphButton("normtog", 0, ID_norm, "Norm", 35, binbox);
			scalebox->burstbutton = scalebox->GraphButton("burstmode", 0, ID_burst, "Burst", 40, filtbox);
			scalebox->selectbutton = scalebox->GraphButton("selectmode", 0, ID_select, "Select", 40, filtbox);
			//scalebox->GraphButton("burstmode", 0, ID_allburst, "All / Select", 70, vbox, 1);
		}	
		scalebox->burstbutton->linkID = ID_select;
		scalebox->selectbutton->linkID = ID_burst;
		vbox->Add(binbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
		vbox->Add(filtbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	}

	if(condex == 3) {
		wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
		if(ostype == Mac) {
			scalebox->ScaleButton(ID_overlay, "Ovl", 43, hbox);
			scalebox->ScaleButton(ID_position, "Pos", 43, hbox);
		}
		else {
			scalebox->ScaleButton(ID_overlay, "Over", 35, hbox);
			hbox->AddSpacer(2);
			scalebox->ScaleButton(ID_position, "Pos", 35, hbox);
		}
		vbox->Add(hbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
		scalebox->overset.Add(ID_overlay, ID_position, 2, 3);
	}

	if(condex == 5) {
		wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
		if(ostype == Mac) {
			scalebox->ScaleButton(ID_overlay2, "Ovl", 43, hbox);
			scalebox->ScaleButton(ID_position2, "Pos", 43, hbox);
		}
		else {
			scalebox->ScaleButton(ID_overlay2, "Over", 35, hbox);
			hbox->AddSpacer(2);
			scalebox->ScaleButton(ID_position2, "Pos", 35, hbox);
		}
		vbox->Add(hbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
		scalebox->overset.Add(ID_overlay2, ID_position2, 4, 5);
	}
}


void MagNetModel::StoreClear()
{
	int i;

	/*
	for(i=0; i<storesize; i++) {
		oxynetdata->water[i] = 0;
	}*/
}


void MagNetModel::RunModel()
{
	if(mainwin->diagnostic) mainwin->SetStatusText("MagNet Model Run");

	ParamStore *netparams = netbox->GetParams();
	numneurons = (*netparams)["numneurons"];

	// Create more neuron objects if requested number is larger than current max
	if(numneurons > modneurons.size()) {
		modneurons.resize((*netparams)["numneurons"]);
		modneurons_max = (*netparams)["numneurons"];
	}

	unsigned long modseed = (*netparams)["modseed"];

	// Set random seed
	if((*netbox->modflags)["seedgen"]) {
		modseed = (unsigned)(time(NULL));
		netbox->paramset.GetCon("modseed")->SetValue(modseed);
	}
	init_mrand(modseed);

	NeuroGen();
	neurodatabox->NeuroData();

    if(!runflag) {
        runflag = true;
        modthread = new MagNetMod(this);
		modthread->diag = false;
        modthread->Create();
        modthread->Run();
    }
    else diagbox->Write("Run command blocked, already running\n\n");
}


void MagNetModel::NeuroGen()
{
	int i, p, numgen, numparams;
	double paramval, paramsdgen;
	double lognormvar;
	wxString *tags;

	//mainwin->diagbox->Write("NeuroGen call\n");
	DiagWrite("NeuroGen call\n");

	ParamStore *netparams = netbox->GetParams();
	ParamStore *neuroparams = spikebox->GetParams();

	//ParamSet *vasoset = vasomod->vasobox->paramset;
	//numgen = vasomod->vasogenbox->numgen;
	//tags = vasomod->vasogenbox->gentags;

	for(i=0; i<numneurons; i++) {
		modneurons[i].type = 0;
		spikebox->GetParams(modneurons[i].spikeparams);
		secbox->GetParams(modneurons[i].secparams);
		signalbox->GetParams(modneurons[i].sigparams);
		dendbox->GetParams(modneurons[i].dendparams);
		synthbox->GetParams(modneurons[i].synthparams);

		if((*netbox->modflags)["netinit"] || !modneurons[i].initflag) {
			// Neuron heterogeneity
			paramsdgen = gaussian(0, 1);
			lognormvar = exp(0 + (*netparams)["synvarsd"] * paramsdgen);
			modneurons[i].synvar = lognormvar;
			// Store initialisation
			modneurons[i].mRNAinit = (*netparams)["mRNAinit"];
			modneurons[i].storeinit = (*netparams)["storeinit"];
		}

		(*modneurons[i].spikeparams)["synvar"] = modneurons[i].synvar;
		(*modneurons[i].synthparams)["mRNAinit"] = modneurons[i].mRNAinit;
		(*modneurons[i].secparams)["Rinit"] = modneurons[i].storeinit;
		modneurons[i].netinit = (*netbox->modflags)["netinit"];
		modneurons[i].storereset = (*netbox->modflags)["storereset"];
		modneurons[i].initflag = true;
	}
}


MagNetModel::~MagNetModel()
{
	delete netdata;
	//delete[] neurons;
	delete currmodneuron;
	delete netdat;
	delete netneuron;
	delete magpop;
	delete neurodata;
}


void MagNetModel::GraphData()
{
	wxString tag, label;
	wxString null = "null";
	int i, selectcode = 100;
	GraphSet *graphset;

	// Data graphs
	//
	// GraphDat(data pointer, xfrom, xto, yfrom, yto, label string, graph type, bin size, colour)
	// ----------------------------------------------------------------------------------
	graphbase->Add(GraphDat(magpop->OxySecretion, 0, 50000, 0, 300, "Oxytocin Secretion", 5, 1, green), "OxySecretion");
	graphbase->Add(GraphDat(magpop->OxyPlasma, 0, 50000, 0, 10000, "Oxytocin Plasma", 5, 1, blue), "OxyPlasma");
	graphbase->Add(GraphDat(&magpop->secLong, 0, 50000, 0, 300, "Secretion 60s", 5, 60, lightblue), "seclong");
	graphbase->Add(GraphDat(&magpop->secHour, 0, 50000, 0, 30, "Secretion 10min", 5, 600, lightblue), "sechour");

	graphbase->Add(GraphDat(&neurodata->secP, 0, 500, 0, 5000, "Secretion P", 5, 1, lightred, 1000/datsample), "oxysecp");
	graphbase->Add(GraphDat(&neurodata->secR, 0, 500, 0, 20000, "Secretion R", 5, 1, lightred, 1000/datsample), "oxysecr");
	graphbase->Add(GraphDat(&neurodata->secX, 0, 500, 0, 30, "Secretion X", 5, 1, lightred, 1000/datsample), "oxysecx");

	graphbase->Add(GraphDat(&neurodata->Ca, 0, 500, 0, 500, "Neuron Ca", 5, 1, lightgreen, 1000/datsample), "neuroCa");

	graphbase->Add(GraphDat(&magpop->OxySecretionNet, 0, 50000, 0, 300, "Net Oxy Secretion", 4, 1, lightgreen), "OxySecretionNet");
	//graphbase->Add(GraphDat(&magpop->secX, 0, 50000, 0, 300, "Net Oxy Secretion", 5, 1, lightgreen), "OxySecretionNet");
	graphbase->Add(GraphDat(&magpop->NetSecretion4s, 0, 50000, 0, 300, "Net Secretion 4s", 5, 4, lightblue), "NetSecretion4s");
	graphbase->Add(GraphDat(&magpop->OxyPlasmaNet, 0, 50000, 0, 10000, "Net Oxy Plasma", 5, 1, lightblue), "OxyPlasmaNet");
	graphbase->Add(GraphDat(&magpop->plasmaLong, 0, 50000, 0, 10000, "Plasma Long", 5, 60, purple), "plasmalong");
	graphbase->Add(GraphDat(&magpop->netsecLong, 0, 50000, 0, 300, "Net Secretion 60s", 5, 60, lightblue), "netseclong");
	graphbase->Add(GraphDat(&magpop->netsecHour, 0, 50000, 0, 300, "Net Secretion 1h", 5, 600, lightblue), "netsechour");

	graphbase->NewSet("Osmotic", "osmo");
	graphbase->GetSet("osmo")->submenu = 1;
	graphbase->Add(GraphDat(&magpop->PlasmaNaConc, 0, 90000, 0, 300, "Plasma Na+ Concentration (1s bins)", 4, 1, lightgreen), "PlasmaNaConc", "osmo");
	graphbase->Add(GraphDat(&magpop->EVFNaConc, 0, 90000, 0, 10000, "EVF Na+ Concentration (1s bins)", 4, 1, lightblue), "EVFNaConc", "osmo");
	graphbase->Add(GraphDat(&magpop->DiffNaGrad, 0, 90000, 0, 300, "Diffusion Na+ Gradient Plasma<->EVF (1s bins)", 4, 1, lightgreen), "DiffNaGrad", "osmo");
	graphbase->Add(GraphDat(&magpop->ICFGrad, 0, 90000, 0, 10000, "ICF<->EVF Gradient (1s bins)", 4, 1, lightblue), "ICFGrad", "osmo");
	graphbase->Add(GraphDat(&magpop->ICFVol, 0, 90000, 0, 300, "ICF Volume (ml) (1s bins)", 4, 1, lightgreen), "ICFVol", "osmo");
	graphbase->Add(GraphDat(&magpop->EVFNaVol, 0, 90000, 0, 10000, "EVF Volume (ml) (1s bins)", 4, 1, lightblue), "EVFNaVol", "osmo");
	graphbase->Add(GraphDat(&magpop->OsmoPress1s, 0, 90000, 0, 10000, "Osmotic Pressure (ml) (1s bins)", 4, 1, lightblue), "OsmoPress1s", "osmo");

	graphbase->Add(GraphDat(&magpop->inputsignal, 0, 50000, 0, 10000, "Input Signal", 5, 1, lightgreen, 10), "inputsignal");
	graphbase->Add(GraphDat(&magpop->netsignal, 0, 50000, 0, 1000, "Net Signal", 5, 1, lightblue), "netsignal");
	graphbase->Add(GraphDat(&magpop->inputLong, 0, 50000, 0, 10000, "Input Long", 5, 60, lightgreen), "inputlong");
	graphbase->Add(GraphDat(&magpop->transLong, 0, 50000, 0, 10000, "Trans Long", 5, 60, lightred), "translong");

	graphbase->Add(GraphDat(&neurodata->pspsig, 0, 50000, 0, 1000, "PSP Signal", 5, 1, lightblue), "pspsig");
	graphbase->Add(GraphDat(&neurodata->V, 0, 50000, 0, 1000, "Rec V", 5, 1, lightgreen), "recV");
	graphbase->Add(GraphDat(&neurodata->syn, 0, 50000, 0, 1000, "Rec Syn", 5, 1, lightblue), "recsyn");
	graphbase->Add(GraphDat(&neurodata->psp, 0, 50000, 0, 1000, "Rec PSP", 5, 1, lightred), "recpsp");
	graphbase->Add(GraphDat(&neurodata->rand, 0, 50000, 0, 1000, "Rec Rand", 5, 1, purple), "recrand");
	//graphbase->Add(GraphDat(&oxyneurodata->inputrate, 0, 50000, 0, 1000, "Input Signal", 5, 1, lightgreen), "inputrate");

	graphbase->Add(GraphDat(&magpop->storesum, 0, 1000, 0, 1000000, "Summed Store", 5, 1000, blue, 1000), "sumstore");
	graphbase->Add(GraphDat(&magpop->storeLong, 0, 1000, 0, 300, "Store Long", 5, 60, lightblue), "storelong");
	graphbase->Add(GraphDat(&magpop->synthstoreLong, 0, 1000, 0, 300, "Synth Store Long", 5, 60, lightred), "synthstorelong");
	graphbase->Add(GraphDat(&magpop->synthrateLong, 0, 1000, 0, 300, "Synth Rate Long", 5, 60, lightred), "synthratelong");
	graphbase->Add(GraphDat(&magpop->storesumLong, 0, 1000, 0, 300, "Summed Store Long", 5, 60, lightblue), "sumstorelong");
	graphbase->Add(GraphDat(&magpop->synthstoresumLong, 0, 1000, 0, 300, "Summed Synth Store", 5, 60, lightred), "sumsynthstorelong");
	graphbase->Add(GraphDat(&magpop->synthratesumLong, 0, 15, 0, 300, "Summed Synth Rate", 5, 60, lightblue), "synthratesumlong");


	graphbase->Add(GraphDat(&neurodata->stimTL, 0, 500, 0, 20000, "Synth TL", 5, 1, lightgreen, 1000/datsample), "stimtl");
	graphbase->Add(GraphDat(&neurodata->stimTS, 0, 500, 0, 30, "Synth TS", 5, 1, lightblue, 1000/datsample), "stimts");
	graphbase->Add(GraphDat(&neurodata->mRNAstore, 0, 500, 0, 30, "mRNA store", 5, 1, lightred, 1000/datsample), "mrnastore");

	graphbase->Add(GraphDat(&magpop->storesumNorm, 0, 1000, 0, 300, "Summed Store Norm", 5, 60, lightblue), "sumstorenorm");


	// Spike data plots

	currmodneuron->PlotSet(graphbase, "Model ", blue, 1, "model");
	netdat->PlotSet(graphbase, "Net ", blue, 1, "net");
	viewcell[0].PlotSet(graphbase, "Cell ", green, 1, "cell0");

	tag = "model";
	graphset = graphbase->NewSet("Model Spikes", "modelspikes");
	graphset->AddFlag("spikeres", 1);
	graphset->AddFlag("rateres", 10);
	graphset->AddFlag("nettog", 100);
	graphset->Add(tag + "rate1s", 0);
	graphset->Add(tag + "spikes1ms", 1);
	graphset->Add(tag + "rate10s", 10);
	graphset->Add(tag + "rate30s", 20);
	graphset->Add(tag + "rate300s", 30);
	graphset->Add("netrate1s", 100);
	graphset->Add("netspikes1ms", 101);
	//graphset->Add("netrate10s", 10);
	if(diagbox) diagbox->textbox->AppendText(graphset->Display());

	graphset = graphbase->NewSet("Model Intervals", "modelintervals");
	graphset->IntervalSet("model");

	/*
	graphset->AddFlag("binrestog1", 1);
	graphset->AddFlag("hazmode1", 10);
	graphset->AddFlag("normtog", 100);
	graphset->AddFlag("nettog", 1000);
	graphset->Add(tag + "hist1ms", 0);
	graphset->Add(tag + "hist5ms", 1);
	graphset->Add(tag + "haz1ms", 10);
	graphset->Add(tag + "haz5ms", 11);
	graphset->Add(tag + "normhist1ms", 100);
	graphset->Add(tag + "normhist5ms", 101);
	graphset->Add(tag + "haz1ms", 110);
	graphset->Add(tag + "haz5ms", 111);
	graphset->Add("nethist1ms", 1000);
	graphset->Add("nethaz1ms", 1010);
	graphset->Add("nethist5ms", 1001);
	graphset->Add("nethaz5ms", 1011);
	*/

	tag = "cell0";
	graphset = graphbase->NewSet("Cell Spikes", "cellspikes");
	graphset->AddFlag("spikeres", 1);
	graphset->AddFlag("rateres", 10);
	graphset->Add(tag + "rate1s", 0);
	graphset->Add(tag + "spikes1ms", 1);
	graphset->Add(tag + "rate100ms", 11);
	graphset->Add(tag + "rate10s", 10);
	graphset->Add(tag + "rate30s", 20);
	graphset->Add(tag + "rate300s", 30);

	graphset = graphbase->NewSet("Cell Intervals", "cellintervals");
	//graphset->IntervalSet("cell0", false, true, selectcode);
	graphset->IntervalSet("cell0");

	//graphbase->GetGraph("cell0spikes1ms")->drawX = 15000;
	//graphbase->GetGraph("cell0spikes1ms")->xaxis = 1;
	//graphbase->GetGraph("cell0spikes1ms")->yaxis = 0;


	tag = "net";
	graphset = graphbase->NewSet("Net Spikes", "netspikes");
	graphset->AddFlag("spikeres", 1);
	graphset->AddFlag("rateres", 10);
	graphset->AddFlag("nettog", 100);
	graphset->Add(tag + "rate1s", 0);
	graphset->Add(tag + "spikes1ms", 1);
	graphset->Add(tag + "rate10s", 10);
	graphset->Add(tag + "rate30s", 20);
	graphset->Add(tag + "rate300s", 30);
	

	//graphbase->Add(GraphDat(&currmodneuron->IoDdata, 0, 100, 0, 100, "IoD Model", 2, 1, lightblue), "iodmod");
	//graphbase->GetGraph("iodmod")->gdatax = &currmodneuron->IoDdataX;
	//graphbase->GetGraph("iodmod")->xcount = 7;  
	//graphbase->GetGraph("iodmod")->synchx = false;

	IoDGraph(&currmodneuron->IoDdata, &currmodneuron->IoDdataX, "IoD Model", "iodmod", lightblue, 0, null);
	IoDGraph(&currmodneuron->burstdata->IoDdata, &currmodneuron->burstdata->IoDdataX, "IoD Burst Model", "iodburstmod", lightblue, 20, null);
	IoDGraph(&currmodneuron->selectdata->IoDdata, &currmodneuron->selectdata->IoDdataX, "IoD Select Model", "iodselectmod", lightgreen, 20, null);
	graphset = graphbase->NewSet("IoD Model", "iodmod");
	graphset->AddFlag("burstmode", 1);
	graphset->AddFlag("selectmode", 10);
	graphset->Add("iodmod", 0);
	graphset->Add("iodburstmod", 1);
	graphset->Add("iodselectmod", 10);
	graphset->Add("iodselectmod", 11);

	graphbase->Add(GraphDat(&nethist, 0, 50, 0, 100, "Net Histogram 1", 1, 1, lightblue), "nethist1");

	IoDGraph(&viewcell[0].IoDdata, &viewcell[0].IoDdataX, "IoD Cell", "iodcell0", lightgreen, 10, null);
	IoDGraph(&viewcell[0].burstdata->IoDdata, &viewcell[0].burstdata->IoDdataX, "IoD Burst Cell", "iodburstcell0", lightblue, 20, null);
	IoDGraph(&viewcell[0].selectdata->IoDdata, &viewcell[0].selectdata->IoDdataX, "IoD Select Cell", "iodselectcell0", lightgreen, 20, null);
	graphset = graphbase->NewSet("IoD Cell", "iodcell0");
	graphset->AddFlag("burstmode", 1);
	graphset->AddFlag("selectmode", 10);
	graphset->Add("iodcell0", 0);
	graphset->Add("iodburstcell0", 1);
	graphset->Add("iodselectcell0", 10);
	graphset->Add("iodselectcell0", 11);

	//graphset->AddFlag("selectmode", 1);
	//graphset->Add("iodcell0", 0);
	//graphset->Add("iodburstcell0", 1);
	//graphset->Add("iodselectcell0", 1);

	/*
	graphbase->Add(GraphDat(oxynetnet->cellpspsig, 0, 1000, 0, 300, "Cell PSP Signal", 4, 1, blue), "cellpspsig");
	graphbase->Add(GraphDat(oxynetnet->cellCa, 0, 1000, 0, 300, "Cell Ca", 4, 1, red, 1), "cellCa");
	graphbase->Add(GraphDat(oxynetnet->cellGnRH, 0, 1000, 0, 300, "Cell GnRH", 4, 1, red, 1), "cellGnRH");
	graphbase->Add(GraphDat(oxynetnet->cellCaLong, 0, 1000, 0, 300, "Cell Ca Long", 4, 1, blue, 1), "cellCaLong");*/

	graphbase->NewSet("Range Data", "rangedata");
	graphbase->GetSet("rangedata")->submenu = 1;

	for(i=0; i<rangecount; i++) {
		tag.Printf("rangedata%d", i);
		label.Printf("Range Data %d", i);
		graphbase->Add(GraphDat(&rangedata[i], 0, 1000, 0, 100, label, 2, 1, lightred + i), tag, "rangedata");
		graphbase->GetGraph(tag)->gdatax = &rangeref;
		graphbase->GetGraph(tag)->xcount = 100;   
		graphbase->GetGraph(tag)->synchx = false;
		graphbase->GetGraph(tag)->scattermode = true;
		graphbase->GetGraph(tag)->scattersize = 5;
	}

	// Data plots
	graphset = graphbase->NewSet("Data Plots", "dataplots");
	graphset->submenu = 1;
	for(i=0; i<10; i++) {
		tag.Printf("datahist%d", i);
		graphbase->Add(GraphDat(&datahist[i], 0, 1000, 0, 100, label.Format("Data Hist%d", i), 1, 1, lightgreen), tag, "dataplots");
		graphbase->GetGraph(tag)->gdatax = &datahistx[i];
		graphbase->GetGraph(tag)->xcount = 5;   
		graphbase->GetGraph(tag)->synchx = false;
		graphbase->GetGraph(tag)->scattermode = true;
		graphbase->GetGraph(tag)->scattersize = 10;
	}

	// Grid data plots
	graphset = graphbase->NewSet("Grid Plots", "gridplots");
	graphset->submenu = 1;
	for(i=0; i<20; i++) {
		tag.Printf("gridplot%d", i);
		graphbase->Add(GraphDat(&gridplot[i], 0, 1000, 0, 100, label.Format("Grid Plot %d", i), 1, 1, lightgreen), tag, "gridplots");
		graphbase->GetGraph(tag)->gdatax = &gridplotx[i];
		graphbase->GetGraph(tag)->xcount = 5;   
		graphbase->GetGraph(tag)->synchx = false;
		graphbase->GetGraph(tag)->scattermode = true;
		graphbase->GetGraph(tag)->scattersize = 10;
		graphbase->GetGraph(tag)->gdataerr = &gridploterr[i];
	}

	gcodes[0] = "modelspikes";
	gcodes[1] = "modelintervals";
	gcodes[2] = "OxySecretion";
	gcodes[3] = "OxySecretionNet";
	gcodes[4] = "OxyPlasmaNet";
	gcodes[5] = "modelspikes";
	gcodes[6] = "modelspikes";
	gcodes[7] = "modelspikes";
	
	gcount = 8;
	gsmode = 1;
}


void MagNetModel::GSwitch(GraphDisp *gpos, ParamStore *gflags)
{
	int i, gdex;
	GraphSet *graphset;
	wxString text;

	// Specify graphs for display

	//mainwin->diagbox->Write("GSwitch call\n");

	if(gsmode == 1) 
		for(i=0; i<gcount; i++) {
			graphset = graphbase->GetSet(gcodes[i]);
			gdex = graphset->GetPlot(i, gflags);
			//if(diagbox) diagbox->textbox->AppendText(text.Format("gpos %d   gcode %s   set %s   plot %d   modesum %d   sdex %d\n", 
			//	i, gcodes[i], graphset->tag, gdex, graphset->modesum, graphset->sdex));
			gpos[i].Front((*graphbase)[gdex]);
			gpos[i].sdex = graphset->sdex;
		}
}


void MagNetModel::DataOutput()
{
	int i, j;
	int row, col;
	int runtime, substeps, modsubsteps;
	wxString text;
	int bincount;

	int binsize = 5;
	int timerange = 1000;
	int gridmax = 500;
	int histcount;
	

	TextGrid *outgrid = gridbox->textgrid[1];
	gridmax = outgrid->GetNumberRows();
	outgrid->CopyUndo();

	ParamStore *calcparams = neurobox->GetParams();
	//viewcell[0].normscale = (*calcparams)["normscale"];
	timerange = (*calcparams)["histrange"];
	if(timerange < 0) timerange = 0;
	histcount = timerange / binsize;
	if(histcount > gridmax) histcount = gridmax;


	// Cell name
	outgrid->SetCell(0, 0, viewcell[0].name);

	// Histogram and Hazard
	outgrid->SetCell(3, 0, "Bin");
	outgrid->SetCell(3, 1, "Hist 5ms");
	outgrid->SetCell(3, 2, "Norm Hist 5ms");
	outgrid->SetCell(3, 3, "Hazard 5ms");
	for(i=0; i<histcount; i++) {
		outgrid->SetCell(i+5, 0, text.Format("%d", i*5));
		outgrid->SetCell(i+5, 1, text.Format("%.0f", viewcell[0].hist5[i]));
		outgrid->SetCell(i+5, 2, text.Format("%.4f", viewcell[0].hist5norm[i]));
		outgrid->SetCell(i+5, 3, text.Format("%.4f", viewcell[0].haz5[i]));
	}

	// Burst Data
	outgrid->SetCell(3, 5, "Bursts");
	outgrid->SetCell(3, 6, text.Format("%d", viewcell[0].burstdata->numbursts));
	outgrid->SetCell(4, 5, "Mean Spikes");
	outgrid->SetCell(4, 6, text.Format("%.2f", viewcell[0].burstdata->meancount));
	outgrid->SetCell(5, 5, "Mean Length");
	outgrid->SetCell(5, 6, text.Format("%.2f", viewcell[0].burstdata->meanlength));
	outgrid->SetCell(6, 5, "Length SD");
	outgrid->SetCell(6, 6, text.Format("%.2f", viewcell[0].burstdata->sdlength));
	outgrid->SetCell(7, 5, "Mean Silence");
	outgrid->SetCell(7, 6, text.Format("%.2f", viewcell[0].burstdata->meansilence));
	outgrid->SetCell(8, 5, "Silence SD");
	outgrid->SetCell(8, 6, text.Format("%.2f", viewcell[0].burstdata->sdsilence));
	outgrid->SetCell(9, 5, "Activity Q");
	outgrid->SetCell(9, 6, text.Format("%.3f", viewcell[0].burstdata->actQ));


	outgrid->SetCell(3, 8, "Rate 30s");
	for(i=0; i<300; i++) {
		outgrid->SetCell(i+5, 8, text.Format("%d", i*30));
		outgrid->SetCell(i+5, 9, text.Format("%.0f", viewcell[0].srate30s[i]));
	}
}
