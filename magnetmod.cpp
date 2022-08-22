/*
*  magnetmod.cpp
*  HypoModel
*
*  Created by Duncan MacGregor
*  University of Edinburgh 2022
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*
*/


#include "magnetmodel.h"

//wxDECLARE_EVENT(wxEVT_COMMAND_MODTHREAD_COMPLETED, wxThreadEvent);


MagNetMod::MagNetMod(MagNetModel *oxynetmod)
	: ModThread(oxynetmod->modbox, oxynetmod->mainwin), neurons(oxynetmod->modneurons) 
{
	mod = oxynetmod;
	spikebox = mod->spikebox;
	synthbox = mod->synthbox;
	netbox = mod->netbox;
	magpop = mod->magpop;

	//neurons = mod->neurons;
	neurodata = mod->neurodata;
	initflag = false;

	//Protocol Parameters
	rampstart = new int[mod->celltypes];
	rampstop = new int[mod->celltypes];
	rampbase = new double[mod->celltypes];
	rampinit = new double[mod->celltypes];
	rampstep = new double[mod->celltypes];
	rampinput = new double[mod->celltypes];
	rampafter = new double[mod->celltypes];
}


void *MagNetMod::Entry()
{
	int i;
	wxString text;

	diagmute = new wxMutex;
	secmute = new wxMutex;
	osmomute = new wxMutex;
    
    //wxCommandEvent endrunevent(wxEVT_COMMAND_TEXT_UPDATED, ID_EndRun);

	Initialise();        // Read in parameters

	// Generate non-independent PSP counts
	if((*netflags)["inputgen"]) {
		mod->netbox->SetStatus("InputGen...\n");
		InputGen();                                   
		mod->netbox->SetStatus("InputGen...OK\n");
	}

	if((*netflags)["realtime"]) {}

	if(prototype == range) RunRange();
	else RunNet();            // Generate and run network and cell threads
	
	//magpop->PopSum();
	SecretionAnalysis();
	//mod->diagbox->Write(text.Format("Pop Sec  Mean %.4f  IoD 1s bin = %.4f\n\n", magpop->secmean, magpop->secIoD));
	//mod->diagbox->Write(text.Format("Pop Sec  Mean %.4f  IoD 4s bin = %.4f\n\n", magpop->secmean_4s, magpop->secIoD_4s));


	// Hetero Pop Analysis    22/2/13

	double synvar, neurate;
	int syndist[1000];
	int ratedist[1000];

	for(i=0; i<1000; i++) {
		syndist[i] = 0;
		ratedist[i] = 0;
	}

	for(i=0; i<numneurons; i++) {
		synvar = (*(neurons[i].spikeparams))["synvar"];
		mainwin->diagbox->Write(text.Format("MagNetMod neuron %d synvar %.4f\n", i, synvar));
		syndist[(int)(synvar*200)/10]++;
		//neurate = neurons[i].ratemean[0];
		neurate = neurons[i].spikecount2 / runtime;
		//ratedist[(int)(neurate*50)/5]++;
		ratedist[(int)(neurate*50)/10]++;   // 0.2 spikes/s bins
	} 

	for(i=0; i<1000; i++) {
		mod->datahistx[0][i] = i * 0.05;
		mod->datahist[0][i] = syndist[i];
		//mod->datahist[0][i] = 100;
		mod->datahistx[1][i] = i * 0.2;
		mod->datahist[1][i] = ratedist[i];
	}

	// Clean Up
	delete diagmute;
	delete secmute;
	delete osmomute;

	 if((*netflags)["inputgen"]) {
		for(int i=0; i<numneurons; i++) {
			delete[] neurons[i].dendinputE;
			delete[] neurons[i].dendinputI;
		}
	}

	delete[] rampstart;
	delete[] rampstop;
	delete[] rampbase;
	delete[] rampinit;
	delete[] rampstep;
	delete[] rampinput;
	delete[] rampafter;

	mod->netbox->SetStatus("");
    
    wxQueueEvent(mod, new wxThreadEvent(EVT_MODTHREAD_COMPLETED));
    //mod->AddPendingEvent(endrunevent);

    //mod->diagbox->Write("Model thread OK\n\n");
    //mainwin->scalebox->GraphUpdate();
    //mod->runflag = false;
    
	return NULL;
}


void MagNetMod::RunRange()
{
	int i, count;
	int rangestart, rangestop, rangestep;
	int rangeindex;
	int inputrate;
	wxString text;
	int secstart, secstop;
	int minstart, minstop;
	double plasmasum, plasmamean;
	double synthsum, synthmean;
	double secsum, secmean;
	double mRNAsum, mRNAmean;
	int startrow = 1;

	ParamStore *protoparams = mod->protobox->GetParams();
	rangestart = (*protoparams)["rangestart"];
	rangestop = (*protoparams)["rangestop"];
	rangestep = (*protoparams)["rangestep"];
	rangeindex = (*protoparams)["rangedata"];

	mainwin->diagbox->Write(text.Format("Runrange %d neurons\n", numneurons));

	count = 0;
	for(inputrate = rangestart; inputrate<=rangestop; inputrate+=rangestep) {
		mod->protobox->currentrange->SetLabel(text.Format("%d", inputrate));    // Display current range parameter
		for(i=0; i<numneurons; i++) (*(neurons[i].spikeparams))["psprate"] = inputrate;  // copy range parameter to neurons
		RunNet();
		// Store plot data
		mod->rangeref[count] = inputrate;
		mod->rangedata[rangeindex][count] = magpop->popfreq;
		// Measure mean plasma
		secstart = 43200;   // 2000;
		secstop = 86400; // 4000;
		minstart = secstart / 60;
		minstop = secstop / 60;

		plasmasum = 0;
		for(i=secstart; i<secstop; i++) plasmasum += magpop->OxyPlasmaNet[i];
		plasmamean = plasmasum / (secstop - secstart);
		mod->rangedata[rangeindex+1][count] = plasmamean;

		secsum = 0;
		for(i=minstart; i<minstop; i++) secsum += magpop->netsecLong[i];
		secmean = secsum / (minstop - minstart);
		mod->rangedata[rangeindex+2][count] = secmean;

		synthsum = 0;
		for(i=minstart; i<minstop; i++) synthsum += magpop->synthratesumLong[i];
		synthmean = synthsum / (minstop - minstart);
		mod->rangedata[rangeindex+3][count] = synthmean;

		mod->gridbox->textgrid[0]->SetCell(count+startrow, 0, text.Format("%.0f", mod->rangeref[count]));
		mod->gridbox->textgrid[0]->SetCell(count+startrow, rangeindex+1, text.Format("%.2f", mod->rangedata[rangeindex][count]));
		mod->gridbox->textgrid[0]->SetCell(count+startrow, rangeindex+2, text.Format("%.2f", mod->rangedata[rangeindex+1][count]));
		mod->gridbox->textgrid[0]->SetCell(count+startrow, rangeindex+3, text.Format("%.2f", mod->rangedata[rangeindex+2][count]));
		mod->gridbox->textgrid[0]->SetCell(count+startrow, rangeindex+4, text.Format("%.2f", mod->rangedata[rangeindex+3][count]));
		count++;
	}

	mod->graphbase->GetGraph(text.Format("rangedata%d", rangeindex))->xcount = count;   
}


void MagNetMod::Initialise()
{
	int i;
	int maxtime = 10000;
	wxString text, tag[10];

	netparams = mod->netbox->GetParams();
		
	netflags = mod->netbox->modflags;

	runtime = int((*netparams)["runtime"]);
	numneurons = int((*netparams)["numneurons"]);
	netrate = int((*netparams)["netrate"]);
	osmorate = int((*netparams)["osmorate"]);
	osmo_hstep = int((*netparams)["osmo_hstep"]);
	buffrate = int((*netparams)["buffrate"]);
	//modseed = (*netparams)["modseed"];
	mod->popscale = (*netparams)["popscale"];

	mod->neurodatabox->neurocount = numneurons;

	secmode = (*netflags)["secmode"];    // run secretion and plasma models if secmode = 1
	plasmamode = (*netflags)["plasmamode"];   

	ParamStore *neuroflags = mod->spikebox->modflags;
	if((*neuroflags)["ipInfusionflag"] || (*neuroflags)["ivInfusionflag"]) osmomode = 1;
	else osmomode = 0;

	/*
	// Set random seed
	if((*netflags)["seedgen"]) {
		modseed = (unsigned)(time(NULL));
		netbox->paramset.GetCon("modseed")->SetValue(modseed);
	}
	init_mrand(modseed);
	*/

	// Initialise osmomod osmotic pressure buffer
	OsmoStore.setsize(maxtime * 1000);
	osmotime = 0;

	// Initialise Population
    magpop->numneurons = numneurons;
	magpop->runtime = runtime;
	magpop->neurons = &neurons;
	//mod->magpop->StoreClear();

	//NeuroGen();     // Copy and generate individual neuron parameters sets

	// Initialise Protocols
	prototype = (*mod->modeflags)["prototype"];

	ParamStore *protoparams = mod->protobox->GetParams();


	if(prototype == ramp) {
		mod->diagbox->Write("prototype ramp\n");
		for(i=0; i<mod->celltypes; i++) {         // currently only one cell type in use, code taken from VMN model
			tag[i].Printf("%d", i);
			rampbase[i] = (*protoparams)["rampbase" + tag[i]];
			rampstart[i] = (*protoparams)["rampstart" + tag[i]];
			rampstop[i] = (*protoparams)["rampstop" + tag[i]];
			rampinit[i] = (*protoparams)["rampinit" + tag[i]];
			if(rampinit[i] < 0) rampinit[i] = rampbase[i];
			rampstep[i] = (*protoparams)["rampstep" + tag[i]];
			rampafter[i] = (*protoparams)["rampafter" + tag[i]];
			if(rampafter[i] < 0) rampafter[i] = rampbase[i] + (rampstop[i] - rampstart[i]) * rampstep[i];
			mainwin->diagbox->Write(text.Format("ramp proto %d  base %.2f  step %.4f  after %.2f\n", i, rampbase[i], rampstep[i], rampafter[i]));
		}
		// Copy proto params to each neuron
		for(i=0; i<numneurons; i++) {
			(*neurons[i].protoparams)["rampbase"] = rampbase[neurons[i].type];
			(*neurons[i].protoparams)["rampstart"] = rampstart[neurons[i].type];
			(*neurons[i].protoparams)["rampstop"] = rampstop[neurons[i].type];
			(*neurons[i].protoparams)["rampinit"] = rampinit[neurons[i].type];
			(*neurons[i].protoparams)["rampstep"] = rampstep[neurons[i].type];
			(*neurons[i].protoparams)["rampafter"] = rampafter[neurons[i].type];
		}
	}

	if(prototype == rampcurve) {
		mod->diagbox->Write("prototype ramp curve\n");
		
		// Copy proto params to each neuron
		for(i=0; i<numneurons; i++) {
			(*neurons[i].protoparams)["rampbase"] = (*protoparams)["rampbase1"];
			(*neurons[i].protoparams)["rampstart"] = (*protoparams)["rampstart1"];
			(*neurons[i].protoparams)["rampstop"] = (*protoparams)["rampstop1"];
			(*neurons[i].protoparams)["rampinit"] = (*protoparams)["rampbase1"];
			(*neurons[i].protoparams)["rampmax"] = (*protoparams)["rampmax1"];
			(*neurons[i].protoparams)["rampgrad"] = (*protoparams)["rampgrad1"];
			(*neurons[i].protoparams)["rampafter"] = (*protoparams)["rampbase1"];
		}
	}

}

/*
void MagNetMod::NeuroGen()
{
	int i, p, numgen, numparams;
	double paramval, paramsdgen;
	double lognormvar;
	wxString *tags;

	//mainwin->diagbox->Write("NeuroGen call\n");
    mod->DiagWrite("NeuroGen call\n");

	ParamStore *netparams = mod->oxynetbox->GetParams();
	ParamStore *neuroparams = mod->spikebox->GetParams();

	//ParamSet *vasoset = vasomod->vasobox->paramset;
	//numgen = vasomod->vasogenbox->numgen;
	//tags = vasomod->vasogenbox->gentags;

	for(i=0; i<numneurons; i++) {
		neurons[i].type = 0;
		mod->spikebox->GetParams(neurons[i].spikeparams);  
		mod->secbox->GetParams(neurons[i].secparams);
		mod->signalbox->GetParams(neurons[i].sigparams);
		mod->dendbox->GetParams(neurons[i].dendparams);
		mod->synthbox->GetParams(neurons[i].synthparams);

		if((*netflags)["netinit"] || !neurons[i].initflag) {
			// Neuron heterogeneity
			paramsdgen = gaussian(0, 1);
			lognormvar = exp(0 + (*netparams)["synvarsd"] * paramsdgen);
			neurons[i].synvar = lognormvar;
			// Store initialisation
			neurons[i].mRNAinit = (*netparams)["mRNAinit"];
			neurons[i].storeinit = (*netparams)["storeinit"];
		}

		(*neurons[i].spikeparams)["synvar"] = neurons[i].synvar;
		(*neurons[i].synthparams)["mRNAinit"] = neurons[i].mRNAinit;
		(*neurons[i].secparams)["Rinit"] = neurons[i].storeinit;
		neurons[i].netinit = (*netflags)["netinit"];
		neurons[i].storereset = (*netflags)["storereset"];
		neurons[i].initflag = true;

		// old heterogeneous bursting parameter code - save for future 13/4/21
		/*vasomod->vasobox->GetParams(params);	    
		for(p=0; p<numgen; p++) {
			paramsdgen = gaussian(0, 1);
			(*params)[tags[p] + "sdgen"] = paramsdgen;
			paramval = (*genparams)[tags[p] + "base"] + (*genparams)[tags[p] + "sd"] * paramsdgen;
			//paramval = gaussian((*genparams)[tags[p] + "base"], (*genparams)[tags[p] + "sd"]);
			if(paramval < 0 && tags[p] != "Vrest") paramval = 0;
			if(tags[p] == "ratioDyno") {
				dynotau = 1/log((double)2)/((*params)["halflifeDyno"]/1000);
				(*params)["kDyno"] = dynotau * paramval; 
			}
			else if(tags[p] == "synvar") {
				lognormvar = exp(0 + (*genparams)[tags[p] + "sd"] * paramsdgen);
				(*params)["synvar"] = lognormvar;
			}
			else (*params)[tags[p]] = paramval;
		}
	}
}
*/


int MagNetMod::InputGen()
{
	int i, n, c, t, inpcell;
	int msteps = 1000;    // number of neuron model steps per s
	double gen;
	double inpfreq;
	int cell;
	int inputcells;
	double neurosyn, synvar;
	double netinput, netIratio;
	double hstep = 1;
	ParamStore *neuroparams;
	wxString text;

	double epspt, ipspt;
	double erate, irate;
	int nepsp, nipsp;
	int numsteps;
	int celltype;

	int maxinputcells = 500;
	int maxconnect = 1000;

	FILE *ofp = NULL, *tofp = NULL;
	wxCommandEvent progevent(wxEVT_COMMAND_TEXT_UPDATED, ID_Progress);

	// std::vector<std::vector<int>> aVector(row_size, std::vector<int>(col_size));

	std::vector<std::vector<short>> Ecellconnect(maxinputcells, std::vector<short>(maxconnect));
	std::vector<std::vector<short>> Icellconnect(maxinputcells, std::vector<short>(maxconnect));
	//short *Icellconnect[500];

	//for(i=0; i<maxinputcells; i++) Ecellconnect[i] = new short[maxconnect];
	//for(i=0; i<maxinputcells; i++) Icellconnect[i] = new short[maxconnect];

	//short Ecellconnect[500][1000];  // Stores the EPSP connection matrix    ([input cell, connection index] = neuron index)  
	//short Icellconnect[500][1000];  // Stores the IPSP connection matrix    max 500 input cells, max 100 neuron connections per input cell
	short Econnect[500];			  // Connection counts for each EPSP input cell 
	short Iconnect[500];			  // Connection counts for each IPSP input cell 
	short inputcheck[500];

	unsigned char *dendinputE[1000];
	unsigned char *dendinputI[1000];

	if(diag) tofp = fopen("inputgen-diag.txt", "w");

	inputcells = (*netparams)["inputcells"];   // Number of input (presynaptic) cells 
	neurosyn = (*netparams)["neurosyn"];       // Number of input (presynaptic) cells connected each neuron in network/population
	netinput = (*netparams)["netinput"];       // EPSP input rate
	netIratio = (*netparams)["netIratio"];     // Input Iratio - IPSP/EPSP ratio

	// inputcells must be greater than or equal to neurosyn
	if(inputcells < neurosyn) {
		mod->diagbox->Write("\nInputGen : Error, inputcells too low for neurosyn\n");
		return 0;
	}

	// Ratio of inputcells to neurosyn determine degree of independence between neurons inputs, 1 gives every neuron the same PSPs 
	// Currently Iratio determines relative IPSP rate rather than number of connected IPSP cells
	// Input rate (netinput) currently fixed but could easily be varied, just have to do rate and Iratio conversion at each step


	// Ramp Protocol
	prototype = (*mod->modeflags)["prototype"];
	double rampstep1ms[2]; 
	rampstep1ms[0] = (double)rampstep[0] / 1000;
	rampstep1ms[1] = (double)rampstep[1] / 1000;
	celltype = 0;

	inpfreq = netinput / neurosyn;
	//inpfreq = 3;
	erate = inpfreq / 1000;
	irate = erate * netIratio;

	//// Set up connections
	//for(i=0; i<inputcells; i++) {
	numsteps = msteps * runtime;

	if(diag) {
		fprintf(tofp, "InputGen %d input cells, %d neurons\n", inputcells, numneurons);
		fprintf(tofp, "Runtime %d  Numsteps %d  Neuron Inputs %.2f\n\n", runtime, numsteps, neurosyn);
		fflush(tofp);
	}

	for(inpcell=0; inpcell<inputcells; inpcell++) {
		Econnect[inpcell] = 0;
		Iconnect[inpcell] = 0;
	}


	// Loop generates random subset of input cells for each neuron, repeated for EPSPs and IPSPs 

	for(n=0; n<numneurons; n++) {
		// Progress
		//if(n%(numneurons/10) == 0) mod->oxynetbox->SetStatus(text.Format("InputGen...%d\%\n", 100*n/numneurons));
		// Input Storage
		neurons[n].dendinputE = new unsigned char[numsteps];
		neurons[n].dendinputI = new unsigned char[numsteps];
		//mod->diagbox->Write(text.Format("neuron %d allocated %d\n"));
		for(t=0; t<numsteps; t++) {
			neurons[n].dendinputE[t] = 0;
			neurons[n].dendinputI[t] = 0;
		}

		// Heterogeneity
		neuroparams = neurons[n].spikeparams;
		synvar = (*neuroparams)["synvar"];

		//fprintf(tofp, "neuron %d storage ok\n", n);
		//fflush(tofp);

		// Connection generation

		for(inpcell=0; inpcell<inputcells; inpcell++) inputcheck[inpcell] = 0;    

		//fprintf(tofp, "neuron %d connections cleared\n", n);
		//fflush(tofp);

		for(c=0; c<(int)(neurosyn * synvar); c++) {
			gen = mrand01();
			cell = floor(gen * inputcells);
			while(inputcheck[cell] > 0) {
				gen = mrand01();
				cell = floor(gen * inputcells);
			}
			inputcheck[cell]++;
			//neurons[n].inputconnect[c] = cell;
			Ecellconnect[cell][Econnect[cell]] = n;
			Econnect[cell]++;
		}

		for(inpcell=0; inpcell<inputcells; inpcell++) inputcheck[inpcell] = 0;

		//fprintf(tofp, "neuron %d connections cleared\n", n);
		//fflush(tofp);

		for(c=0; c<(int)(neurosyn * synvar); c++) {
			gen = mrand01();
			cell = floor(gen * inputcells);
			while(inputcheck[cell] > 0) {
				gen = mrand01();
				cell = floor(gen * inputcells);
			}
			inputcheck[cell]++;
			//neurons[n].inputconnect[c] = cell;
			Icellconnect[cell][Iconnect[cell]] = n;
			Iconnect[cell]++;
		}

		//fprintf(tofp, "neuron %d connections generated\n", n);
		//fflush(tofp);
	}

	if(diag) fprintf(tofp, "\n");

	mod->netbox->SetStatus("InputGen...connect OK...\n");

	// Input generation

	// Loop generates random epsps for each input cell
	// - could reverse loop order to put cell loop inside time step loop for more efficient variable erate
	// - but current order better for potential multi-threading

	for(inpcell=0; inpcell<inputcells; inpcell++) {
		epspt = -log(1 - mrand01()) / erate;   // initial epspt value, could just be 0
		ipspt = -log(1 - mrand01()) / irate;

		mod->netbox->SetStatus(text.Format("InputGen...connect OK...cell %d\n", inpcell));

		// input array pointer copies for fast access
		for(c=0; c<Econnect[inpcell]; c++) dendinputE[c] = neurons[Ecellconnect[inpcell][c]].dendinputE;
		for(c=0; c<Iconnect[inpcell]; c++) dendinputI[c] = neurons[Icellconnect[inpcell][c]].dendinputI;

		for(t=0; t<numsteps; t++) {

			// Variable Input Signal
			if(prototype == ramp) {
				//mod->diagbox->Write("set ramp\n"); 
				if(t < rampstart[celltype]*1000) rampinput[celltype] = rampbase[celltype];
				if(t >= rampstart[celltype]*1000 && t < rampstop[celltype]*1000) 
					rampinput[celltype] = rampinit[celltype] + (t - rampstart[celltype]*1000) * rampstep1ms[celltype] * hstep;
				//if(t >= rampstop[celltype]*1000) rampinput[celltype] = rampbase[celltype];
				if(t >= rampstop[celltype]*1000) rampinput[celltype] = rampafter[celltype];
				netinput = rampinput[celltype];
				if(netinput < 0) netinput = 0;

				inpfreq = netinput / neurosyn;
				erate = inpfreq / 1000;
				irate = erate * netIratio;
			}

			if(!inpcell && t%1000 == 0) magpop->netsignal[t/1000] = netinput;   // Record net signal at 1s resolution

			nepsp = 0;
			while(epspt < hstep) {
				nepsp++;
				epspt = -log(1 - mrand01()) / erate + epspt;
			}
			epspt = epspt - hstep;

			// At each step add EPSPs to EPSP counts for all connected neurons
			if(nepsp > 0) {
				// add input
				for(c=0; c<Econnect[inpcell]; c++) 
					//neurons[Ecellconnect[inpcell][c]].dendinputE[t] += nepsp;	
					dendinputE[c][t] += nepsp;
			}

			nipsp = 0;
			while(ipspt < hstep) {
				nipsp++;
				ipspt = -log(1 - mrand01()) / irate + ipspt;
			}
			ipspt = ipspt - hstep;
			if(nipsp > 0) {
				// add input
				for(c=0; c<Iconnect[inpcell]; c++)
					//neurons[Icellconnect[inpcell][c]].dendinputI[t] += nipsp;
					dendinputI[c][t] += nipsp;
			}
		}
	}

	/*for(inpcell=0; inpcell<inputcells; inpcell++) {
		ipspt = -log(1 - mrand01()) / irate;
		for(t=0; t<numsteps; t++) {
			nipsp = 0;
			while(ipspt < hstep) {
				nipsp++;
				ipspt = -log(1 - mrand01()) / irate + ipspt;
			}
			ipspt = ipspt - hstep;
			if(nipsp > 0) {
				// add input
				for(c=0; c<Iconnect[inpcell]; c++)
					neurons[Icellconnect[inpcell][c]].dendinputI[t] += nipsp;
			}
		}
	}*/

	// Diagnostic Output

	if(diag) {
		ofp = fopen("inputnet.txt", "w");
		
		fprintf(ofp, "Input Network\n\n");
		fprintf(ofp, "%d input cells\n", inputcells);

		for(i=0; i<10; i++) {
			fprintf(ofp, "input cell %d  connections %d\n", i, Econnect[i]);
			fprintf(ofp, "neuron ");
			for(c = 0; c<Econnect[i]; c++) fprintf(ofp, "%d ", Ecellconnect[i][c]);
			fprintf(ofp, "\n");
		}
		fprintf(ofp, "\n"); 

		int epspsum = 0, ipspsum = 0;

		for(t=0; t<1000; t++) {
			if(neurons[0].dendinputE[t] > 0 || neurons[0].dendinputI[t] > 0) 
				fprintf(ofp, "neuron0  time %dms  epsp %d  ipsp %d\n", t, neurons[0].dendinputE[t], neurons[0].dendinputI[t]);
			epspsum += neurons[0].dendinputE[t];
			ipspsum += neurons[0].dendinputI[t];
		}
		fprintf(ofp, "epsp total %d  ipsp total %d\n", epspsum, ipspsum);

		fclose(ofp);
	}
	if(diag) fclose(tofp);

	return 1;
}


void MagNetMod::RunNet()
{
	int i;
	int step;
	wxString text;
	int maxtime = magpop->maxtime;
	clock_t timestart, timerun;

	mod->DiagWrite(text.Format("\nRunNet %d neurons\n\n", numneurons));

	netsecX = 0;
	tPlasma = 0;
	tEVF = 0;

	// Initialise buffered secretion summation store
	for(i=0; i<maxtime*1000; i++) mod->magpop->secX[i] = 0;
	for(i=0; i<maxtime; i++) mod->magpop->secXcount[i] = 0;
	mod->magpop->secXtime = -1;

	// Generate and run neuron threads
	// Every thread is an instance of the class MagNeuroMod that runs the single neuron code 
	// Every thread needs to be created, run, and deleted after it has finished

	// Create Threads
	for(i=0; i<numneurons; i++) {
		//mod->diagbox->Write(text.Format("Init cell %d\n", i));
		neurothread[i] = new MagNeuroMod(i, &neurons[i], this); 
		neurothread[i]->Create();
	}
	//if(osmomode) osmothread = new OxyOsmoMod(this);
	if(plasmamode) plasmathread = new MagPlasmaMod(this);

	timestart = clock();

	// Run Threads
	for(i=0; i<numneurons; i++) neurothread[i]->Run(); 
	//if(osmomode) osmothread->Run();
	if(plasmamode) plasmathread->Run();

	// Wait for Thread Completion
	for(i=0; i<numneurons; i++) {
		neurothread[i]->Wait(); 
		//mod->diagbox->Write(text.Format("Cell %d OK\n", i));
	}
	//if(osmomode) osmothread->Wait();
	if(plasmamode) plasmathread->Wait();

	timerun = clock() - timestart;
	if(mainwin->diagbox) mod->DiagWrite(text.Format("runtime %d clicks (%f seconds)\n", timerun,((double)timerun)/CLOCKS_PER_SEC));

	// Clean up threads
	for(i=0; i<numneurons; i++) delete neurothread[i]; 
	//if(osmomode) delete osmothread;
	if(plasmamode) delete plasmathread;

	mod->DiagWrite(text.Format("\nRunNet OK\n\n"));

	mod->neurodatabox->NeuroData();


	//
	// Population Analysis
	//

	//double popfreq = 0;  // Global FR

	magpop->PopSum();

	/*for(i=0; i<numneurons; i++) {
		neurons[i].ratecalc();
		popfreq += neurons[i].freq;
	}*/

	//popfreq = popfreq / numneurons;

	

	if((*netflags)["netanalysis"]) {

		// 'netspikes' is a SpikeDat used to store population summed/mean neuron analysis 
		// 'netneuron' is a SpikeDat used for temporary individual neuron analysis
	
		// Population spike rate, histogram, and hazard in different binwidths 

		// Reset population counts
		for(i=0; i<magpop->maxtime; i++) mod->netdat->srate1s[i] = 0;
		for(i=0; i<magpop->maxtime/10; i++) mod->netdat->srate10s[i] = 0;
		for(i=0; i<magpop->maxtime/30; i++) mod->netdat->srate30s[i] = 0;
		for(i=0; i<magpop->maxtime/300; i++) mod->netdat->srate300s[i] = 0;
		for(i=0; i<magpop->maxtime/600; i++) mod->netdat->srate600s[i] = 0;

		for(i=0; i<1000000; i++) mod->netdat->srate1[i] = 0;   // 1ms bins for individual spikes

		for(i=0; i<10000; i++) {
			mod->netdat->hist1[i] = 0;
			mod->netdat->hist5[i] = 0;
			mod->netdat->haz1[i] = 0;
			mod->netdat->haz5[i] = 0;
		}

		magpop->popfreq = 0;

		// Analyse and sum each neuron
		for(i=0; i<numneurons; i++) {
			mod->netneuron->neurocalc(&(neurons[i]));
			magpop->popfreq += mod->netneuron->freq;

			for(step=0; step<magpop->maxtime; step++) mod->netdat->srate1s[step] += mod->netneuron->srate1s[step];  // 1s bins
			for(step=0; step<magpop->maxtime/10; step++) mod->netdat->srate10s[step] += mod->netneuron->srate10s[step];  // 10s bins
			for(step=0; step<magpop->maxtime/30; step++) mod->netdat->srate30s[step] += mod->netneuron->srate30s[step];  // 30s bins
			for(step=0; step<magpop->maxtime/300; step++) mod->netdat->srate300s[step] += mod->netneuron->srate300s[step];  // 300s bins
			for(step=0; step<magpop->maxtime/600; step++) mod->netdat->srate600s[step] += mod->netneuron->srate600s[step];  // 600s bins

			for(step=0; step<1000000; step++) mod->netdat->srate1[step] += mod->netneuron->srate1[step];  // 1 ms bins

			for(step=0; step<10000; step++) {
				mod->netdat->hist1[step] += mod->netneuron->hist1[step];
				mod->netdat->hist5[step] += mod->netneuron->hist5[step];
				mod->netdat->haz1[step] += mod->netneuron->haz1[step];
				mod->netdat->haz5[step] += mod->netneuron->haz5[step];
			}
		}

		// Convert sums to means
		for(step=0; step<magpop->maxtime; step++) mod->netdat->srate1s[step] = mod->netdat->srate1s[step] / numneurons;  // 1s bins
		for(step=0; step<magpop->maxtime/10; step++) mod->netdat->srate10s[step] = mod->netdat->srate10s[step] / numneurons;  // 10s bins
		for(step=0; step<magpop->maxtime/30; step++) mod->netdat->srate30s[step] = mod->netdat->srate30s[step] / numneurons;  // 30s bins
		for(step=0; step<magpop->maxtime/300; step++) mod->netdat->srate300s[step] = mod->netdat->srate300s[step] / numneurons;  // 300s bins
		for(step=0; step<magpop->maxtime/600; step++) mod->netdat->srate600s[step] = mod->netdat->srate600s[step] / numneurons;  // 600s bins

		magpop->popfreq = magpop->popfreq / numneurons;
	}

	mod->DiagWrite(text.Format("\n%d neurons   pop freq %.4f\n", numneurons, magpop->popfreq));
}


void MagNetMod::Export2file(int steps, wxString filename, datdouble vector2print)
{
	float tempvalue;
	wxTextFile stfile; //Secretion rate file
	wxString stfilename = "E:/Data/Test/Average";
	stfilename = stfilename.Append(wxString::Format(wxT("%s"), filename));
	stfilename = stfilename.Append(wxString::Format(wxT("%d"), numneurons));
	stfilename = stfilename.Append("neurones.dat");
	if (!stfile.Create(stfilename)) stfile.Open(stfilename);
	stfile.Clear();
	tempvalue = 0;
	for (double step=0; step<steps; step++) 
	{
		tempvalue = (float)vector2print[step];// / (float)numneurons;  // averaging to a floating point number  
		stfile.AddLine(wxString::Format(wxT("%.2f"), tempvalue ));
	}
	stfile.Write(); // necessary for saving changes! 
	stfile.Close();  // Closing the file					
}


void MagNetMod::SecretionAnalysis()
{
	int i, datacount;
	double mean, variance;

	// secretion rate IoD - 1s bin
	mean = 0;
	variance = 0;
	for(i=0; i<runtime; i++) mean = mean + mod->magpop->OxySecretionNet[i];	   // mean
	mean = mean / runtime;
	for(i=0; i<runtime; i++) variance += (mean - mod->magpop->OxySecretionNet[i]) * (mean - mod->magpop->OxySecretionNet[i]);	  // variance
	variance = variance / runtime;
	mod->magpop->secIoD = variance / mean;
	mod->magpop->secmean = mean;

	// secretion rate IoD - 4s bin
	mean = 0;
	variance = 0;
	datacount = runtime / 4;
	for(i=0; i<datacount; i++) mean = mean + magpop->NetSecretion4s[i];	   // mean
	mean = mean / datacount;
	for(i=0; i<datacount; i++) variance += (mean - magpop->NetSecretion4s[i]) * (mean - magpop->NetSecretion4s[i]);	  // variance
	variance = variance / datacount;
	magpop->secIoD_4s = variance / mean;
	magpop->secmean_4s = mean;
}


