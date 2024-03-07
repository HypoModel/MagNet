/*
*  magneuromod.cpp
*  HypoModel
* 
*  Created by Duncan MacGregor
*  University of Edinburgh 2022
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*  Created June 2017
*  Modified November 2017
*
*/


#include "magnetmodel.h"
#include <random>
#include <math.h>



inline double vox_tanh( const double x )
{
	const double ax = fabs( x );
	const double x2 = x * x;
	const double z = x * ( 1.0 + ax +
		( 1.05622909486427 + 0.215166815390934 * x2 * ax ) * x2 );

	return( z / ( 1.02718982441289 + fabs( z )));
}


MagNeuroMod::MagNeuroMod(int index, MagNeuron *oxyneuron, MagNetMod *oxynetmod)
	: wxThread(wxTHREAD_JOINABLE)
{
	wxString text;

	neuron = oxyneuron; // individual neuron pointer
	netmod = oxynetmod;  

	mod = netmod->mod;
	magpop = mod->magpop;
	neurodex = index;  // setting the index for the current neuron run
	diagbox = netmod->mod->diagbox;
	netbox = netmod->netbox;
	neurorecord = mod->neurodata;

	maxtime = magpop->maxtime;
	maxtimeLong = magpop->maxtimeLong;

	//net->mod->diagbox->Write(text.Format("Cell %d initialising\n", celldex));

	ParamStore *spikeparams = neuron->spikeparams;
	ParamStore *secparams = neuron->secparams;
	ParamStore *sigparams = neuron->sigparams;
	ParamStore *dendparams = neuron->dendparams;
	ParamStore *synthparams = neuron->synthparams;
	ParamStore *protoparams = neuron->protoparams;

	ParamStore *neuroflags = netmod->spikebox->modflags;
	ParamStore *synthflags = netmod->synthbox->modflags;
	ParamStore *netparams = netbox->GetParams();

	modsteps = oxynetmod->runtime * 1000;
	hstep = (*spikeparams)["hstep"];
	shstep = hstep / 1000;
	syn_hstep = shstep / 3600;
	netrate = oxynetmod->netrate;
	osmorate = oxynetmod->osmorate;
	buffrate = oxynetmod->buffrate;

	modseed = (*netparams)["modseed"];
	disprate = (*netparams)["disprate"];

	// Flags
	ipInfusionflag = (*neuroflags)["ipInfusionflag"];
	ivInfusionflag = (*neuroflags)["ivInfusionflag"];
	if(ipInfusionflag || ivInfusionflag) osmomode = 1;
	else osmomode = 0;
	epspsynchflag = (*neuroflags)["epspsynchflag"];
	modmode = (*neuroflags)["modmode"];
	AHP2mode = (*neuroflags)["AHP2mode"];   // 0 for basic, 1 for Ca threshold AHP (PLoS 2012)
	dynostoreflag = (*neuroflags)["dynostoreflag"];

	// Input Parameters
	noimean = (*sigparams)["noimean"];
	noitau = (*sigparams)["noitau"];
	noiamp = (*sigparams)["noiamp"];
	sigIratio = (*sigparams)["sigIratio"];
	signalmode = (*netmod->mod->signalbox->modflags)["noiseflag"];

	// Spiking
	Vthresh = (*spikeparams)["Vthresh"];
	Vrest = (*spikeparams)["Vrest"];
	pspmag = (*spikeparams)["pspmag"];
	psprate = (*spikeparams)["psprate"];
	iratio = (*spikeparams)["iratio"];
	halflifeMem = (*spikeparams)["halflifeMem"];
	kHAP = (*spikeparams)["kHAP"];
	halflifeHAP = (*spikeparams)["halflifeHAP"];
	kDAP = (*spikeparams)["kDAP"];
	halflifeDAP = (*spikeparams)["halflifeDAP"];
	kAHP = (*spikeparams)["kAHP"];
	halflifeAHP = (*spikeparams)["halflifeAHP"];
	synvar = (*spikeparams)["synvar"];

	if(modmode == 1) {
		kAHP2 = (*spikeparams)["kAHP2"];
		halflifeAHP2 = (*spikeparams)["halflifeAHP2"];
		aAHP2 = (*spikeparams)["aAHP2"];
		kDyno = (*spikeparams)["kDyno"];
		halflifeDyno = (*spikeparams)["halflifeDyno"];
		kCa = (*spikeparams)["kCa"];
		halflifeCa = (*spikeparams)["halflifeCa"];
		Ca_rest = (*spikeparams)["Ca_rest"];
		ka = (*spikeparams)["ka"];
		gKL = (*spikeparams)["gKL"];
		gOsmo = (*spikeparams)["gOsmo"];
	}
	else {
		kAHP2 = 0;
		halflifeAHP2 = 1000;
		aAHP2 = 0;
		kDyno = 0;
		halflifeDyno = 1000;
		kCa = 0;
		halflifeCa = 1000;
		Ca_rest = 0;
		ka = 0;
		gKL = 0;
		gOsmo = 0;
	}

	pspmag2 = (*spikeparams)["pspmag2"];
	psprate2 = (*spikeparams)["psprate2"];
	halflifePSP2 = (*spikeparams)["halflifePSP2"];

	
	/*// NMDA EPSP
	if(modmode == 0) {
		pspmag2 = (*spikeparams)["pspmag2"];
		psprate2 = (*spikeparams)["psprate2"];
		halflifePSP2 = (*spikeparams)["halflifePSP2"];
	}
	else {
		pspmag2 =  0;
		psprate2 = 0;
		halflifePSP2 = 1000;
	}*/


	// Secretion
	kB = (*secparams)["kB"];
	halflifeB = (*secparams)["halflifeB"];
	Bbase = (*secparams)["Bbase"];
	kC = (*secparams)["kC"];
	halflifeC = (*secparams)["halflifeC"];
	kE = (*secparams)["kE"];
	halflifeE = (*secparams)["halflifeE"];
	Cthresh = (*secparams)["Cth"];
	Cgradient = (*secparams)["Cgradient"];
	Ethresh = (*secparams)["Eth"];
	Egradient = (*secparams)["Egradient"];
	alpha = (*secparams)["alpha"] / 1000;
	beta = (*secparams)["beta"] / 1000;
	Rmax = (*secparams)["Rmax"];
	Rinit = (*secparams)["Rinit"];
	Pmax = (*secparams)["Pmax"];
	secExp = (*secparams)["secExp"];
	plasma_hstep = (*secparams)["plasma_hstep"];

	// Synthesis
	//vsynrate = (*synthparams)["vsynrate"];  
	vsynrate = 0;       // no longer in use 10/8/21
	// vtrans = (*synthparams)["vtrans"]; 
	vtrans = 0;  // no longer in use 26/6/22
	mRNAmax = (*synthparams)["mRNAmax"];
	mRNAinit = (*synthparams)["mRNAinit"];
	mRNAhalflife = (*synthparams)["mRNAhalflife"];
	synthdel = (*synthparams)["synthdel"]; 
	halflifeTS = (*synthparams)["halflifeTS"];
	kTS = (*synthparams)["kTS"];
	maxTL = (*synthparams)["maxTL"];
	kTL = (*synthparams)["kTL"];
	halflifeTL = (*synthparams)["halflifeTL"];
	basalTL = (*synthparams)["basalTL"];
	rateSR = (*synthparams)["rateSR"];
	synscale = (*synthparams)["synscale"];
	transflag = (*synthflags)["transflag"];
	synthmode = (*synthflags)["synthmode"];
	transmode = (*synthflags)["transmode"];
	scalemode = (*synthflags)["scalemode"];
	polymode = (*synthflags)["polymode"];
	decaymode = (*synthflags)["decaymode"];

	// Dynamic Dynorphin Parameters
	kstoreDyno = (*dendparams)["kstoreDyno"];
	halflifestoreDyno = (*dendparams)["halflifestoreDyno"];
	spikeDyno = (*dendparams)["spikeDyno"];
	tauDynoup = (*dendparams)["tauDynoup"];
	kdendCa = (*dendparams)["kdendCa"];
	halflifedendCa = (*dendparams)["halflifedendCa"];

	// Osmotic Pressure
	BasalNaConc = (*spikeparams)["BasalNaConc"];  // in mOsmoles/l

	// Diffusion and Clearance
	halflifeClear = (*secparams)["ClearHL"];
	halflifeDiff = (*secparams)["DiffHL"];
	PlasmaVol = (*secparams)["VolPlasma"];
	EVFVol = (*secparams)["VolEVF"];

	// Protocol
	prototype = (*mod->modeflags)["prototype"];
	rampbase = (*protoparams)["rampbase"];
	rampstart = (*protoparams)["rampstart"] * 1000;
	rampstop = (*protoparams)["rampstop"] * 1000;
	rampinit = (*protoparams)["rampinit"];
	rampstep = (*protoparams)["rampstep"] / 1000;
	rampafter = (*protoparams)["rampafter"];
	rampmax = (*protoparams)["rampmax"];
	rampgrad = (*protoparams)["rampgrad"] / 1000000;       // scaled for vaso synth
}


void *MagNeuroMod::Entry()
{
	wxString text;

	/*net->diagmute->Lock();
	net->mod->diagbox->Write(text.Format("Cell %d running\n", celldex));
	net->diagmute->Unlock();*/

	neuromod();

	/*net->diagmute->Lock();
	net->mod->diagbox->Write(text.Format("Cell %d finished\n", celldex));
	net->diagmute->Unlock();*/

	return NULL;
}


// Neural model code including integrated spiking, secretion, and synthesis models

void MagNeuroMod::neuromod() 
{
	int i, step;
	int runtime, runtime100;
	wxString text;
	unsigned long seed; 
	double erand, irand;
	int maxspikes;
	bool flagError = false;
	//sfmt_t sfmt;   // new SFMT random number generator  July 2020

	double epspt, ipspt;
	double pspRatio;
	int nepsp, nipsp;
	double epspt1, ipspt1;
	int nepsp1, nipsp1;

	double inputPSP, inputPSP1;
	double ttime, neurotime;
	bool monitor = true;
	bool countflag = false;

	double epsprate, totalepsprate, epspmag;
	double ipsprate, totalipsprate, ipspmag;
	double tauMem;
	double tauHAP;
	double tauDAP;
	double tauAHP, tauAHP2;
	double tauCa, tauDyno;
	double taudendCa;
	double absref;
	double IKL, KLact;

	double tauB, tauE, tauC;
	double tauClear, tauDiff;

	// NMDA synapse EPSPs - new February 2020
	double epspt2;        
	int nepsp2;
	double epsprate2, epspmag2;
	double inputPSP2;
	double tauPSP2;

	double inputOsmo;

	// Variables
	double pspsig;
	double V;
	double tCa, tdendCa;
	double tHAP;
	double tDAP;
	double tAHP, tAHP2;
	double tDyno, storeDyno;

	double tB;
	double tE, tC;
	double tR, tP;
	double tOxyPlasma, tOxyEVF;
	double CaEnt, Cinh, Einh;
	double secX, secBinX, DiffRate;
	double oldsecX;   // used for updating summed population secretion
	double netsecX; 
	double EKpow, CKpow;
	double Ethpow, Cthpow;

	double secRate1s, plasmaRate1s;
	double secRate60s, secRate600s;       // new 24/2/21 for vaso synth figures
	double netsecRate1s, netplasmaRate1s;

	double OsmoPress;  // Instantaneous osmotic pressure
	double IrOsmoPress;  // Instantaneous PSP change due to the osmotic pressure
	double OsmoSetPoint, OsmoShift;  // Set point for homeostatic osmolality. Set at 302,5 m-osmole/kg

	int buffdex;
	double *secXbuffer = new double[buffrate];
	double *secXpop = magpop->secX.data.data();

	double synsig, noisig;
	double epsprate1, ipsprate1;
	int inputgen;
	double rampinput;

	clock_t timestart, timerun;
	int recstart, recstop;        // time points for diagnostic recording
	int recneuron;

	int synthdex;
	int synthrecrate = 1000 * 60;

	maxspikes = neuron->maxspikes;

	int datsample = netmod->mod->datsample;
	//if(celldex == netmod->currentcell) countflag = true; 
	if(neurodex == 0) {
		countflag = true; 
		//net->mod->diagbox->Write(text.Format("cellmod %d\n", celldex));
	}

	wxCommandEvent progevent(wxEVT_COMMAND_TEXT_UPDATED, ID_Progress);
	wxCommandEvent plotevent(wxEVT_COMMAND_TEXT_UPDATED, ID_Display);

	runtime = netmod->runtime * 1000;
	runtime100 = runtime / 100;

	inputgen = (*netmod->netflags)["inputgen"];

	//FILE *tofp = fopen("oxynetneuromod.txt", "w");


	// Parameters
	// If there is NaCl infusion, we interpret that all the PSP are due to the osmotic pressure
	if(osmomode) epsprate = 0;
	else epsprate = psprate / 1000;

	pspRatio = iratio;
	ipsprate = epsprate * pspRatio;
	epspmag = pspmag;
	ipspmag = pspmag;
	epspmag2 = pspmag2;
	epsprate2 = psprate2 / 1000;
	absref = 2;

	// Time Constants - conversion from half-life
	// Spiking 
	tauMem = log((double)2) / halflifeMem;
	tauHAP = log((double)2) / halflifeHAP;
	tauDAP = log((double)2) / halflifeDAP;
	tauAHP = log((double)2) / halflifeAHP;
	tauAHP2 = log((double)2) / halflifeAHP2;

	// Dendritic
	tauCa = log((double)2) / halflifeCa;
	tauDyno = log((double)2) / halflifeDyno;
	taudendCa = log((double)2) / halflifedendCa;

	// NMDA PSP
	tauPSP2 = log((double)2) / halflifePSP2;

	// Secretion
	tauB = log((double)2) / halflifeB;
	tauC = log((double)2) / halflifeC;
	tauE = log((double)2) / halflifeE;

	// Plasma
	tauClear = log((double)2) / (halflifeClear * 1000);
	tauDiff = log((double)2) / (halflifeDiff * 1000);

	// Synthesis
	tauTS = log((double)2) / halflifeTS;
	tauTL = log((double)2) / halflifeTL;
	mRNAtau = log((double)2) / mRNAhalflife;

	double *synthrec = new double[35000];  // record synthrate for delay recall, minute sampled capacity for 24 days


	// initialise random number generator
	seed = modseed + neurodex;
	//seed = 1568637350;
	//para_init_mrand(neurodex, seed);
	//sfmt_init_gen_rand(&sfmt, seed);

	std::mt19937 randgen(seed);
	std::uniform_real_distribution<float> unif01(0, 1);



	// random number test
	
	//for(i=0; i<10; i++) {
	//	erand = sfmt_genrand_real2(&sfmt);
	//	mod->diagbox->Write(text.Format("SFMT random %.4f\n", erand));
	//}

	//std::mt19937 randmt (seed);
	//std::uniform_real_distribution<double> dis (0.0, 1.0);
	//double randomRealBetweenZeroAndOne = dis(generator);

	//seed = (unsigned)(time(NULL));
	//para_init_mrand(neurodex, seed + neurodex);  // is it starting randomly from one of the neurones?

	// Initialise
	totalepsprate = 0;
	totalipsprate = 0;
	epspt = 0;
	ipspt = 0;
	epspt1 = 0;
	ipspt1 = 0;
	pspsig = 0;
	ttime = 0;
	neurotime = 0;
	tHAP = 0;
	tDAP = 0;
	tAHP = 0;
	tAHP2 = 0;
	tCa = Ca_rest;
	tdendCa = 0;
	tDyno = 0;
	storeDyno = 0.6;
	V = Vrest;

	// NMDA PSP
	epspt2 = 0;
	inputPSP2 = 0;

	// Osmotic pressure
	OsmoSetPoint = BasalNaConc * 2;
	OsmoPress = OsmoSetPoint;
	netmod->OsmoPress = OsmoSetPoint;
	if(osmomode) IrOsmoPress = (26 * (OsmoPress - 303)) / 1000; // differential PSP due to the hyperosmotic injection
	else IrOsmoPress = 0;

	/*if (OsmoTimeIv*1000 < modsteps) FlagOsmoIv = true;
	else FlagOsmoIv = false;*/

	// Secretion model
	tR = Rinit;  // Reserve Pool
	tP = Pmax;  // Releasable Pool 
	tOxyPlasma = 0;
	tOxyEVF = 0;
	tB = 0;  // Broadening
	tE = 0;  // Fast Ca2+
	tC = 0.03; // Slow Ca2+

	Ethpow = Ethresh * Ethresh * Ethresh * Ethresh * Ethresh;  // precalculate instead of each loop
	Cthpow = Cthresh * Cthresh * Cthresh;

	// Synthesis
	stimTS = 0;
	stimTL = 0;
	mRNAstore = mRNAinit;
	synthrec[0] = 0;

	// Population secretion and plasma
	secBinX = 0;
	secRate1s = 0;
	secRate60s = 0;
	secRate600s = 0;
	plasmaRate1s = 0;
	netsecRate1s = 0;     // For neuron 0 recording population secretion and plasma in 1s window
	netplasmaRate1s = 0;
	secX = 0;
	oldsecX = 0;
	buffdex = 0;

	neuron->spikecount = 0;
	neuron->spikecount2 = 0;

	noisig = noimean;

	for(double i=0; i<(modsteps/1000); i++) {
		neuron->Secretion[i] = 0;
		//neuron->OxyPlasma[i] = 0;		
	}

	// Record Initial Values
	neuron->storeLong[0] = tR;
	neuron->transLong[0] = 0;
	neuron->synthstoreLong[0] = mRNAstore;
	neuron->synthrateLong[0] = rateSR * (stimTL + basalTL) * synscale * mRNAstore * 3600;
	neurorecord->stimTL[0] = stimTL;
	neurorecord->stimTS[0] = stimTS;
	neurorecord->mRNAstore[0] = mRNAstore;
	neurorecord->Ca[0] = Ca_rest;
	magpop->inputsignal[0] = psprate;
	if(prototype == ramp || prototype == rampcurve) magpop->inputLong[0] = rampbase;
	else magpop->inputLong[0] = psprate;

	// Diagnostic
	if(netmod->diag && neurodex == 0) {
		TextFile diagfile;
		diagfile.New("oxynetcelldiag.txt");
		diagfile.WriteLine(text.Format("oxynet cell %d running for %d steps\n", neurodex, modsteps));
		diagfile.WriteLine(text.Format("Vrest %.2f  Vthresh %.2f\n", Vrest, Vthresh));
		diagfile.WriteLine(text.Format("tauMem %.4f  tauHAP %.4f  tauAHP %.4f\n", tauMem, tauHAP, tauAHP));
		diagfile.WriteLine(text.Format("kCa %.4f  tauCa %.4f  kDyno %.4f  tauDyno %.4f\n", kCa, tauCa, kDyno, tauDyno));
		diagfile.Close();
	}

	//fprintf(tofp, "seed %lu\n", seed);

	//if(prototype == ramp) mod->diagbox->Write(text.Format("ramp base %.2f step %.4f\n", rampbase, rampstep));

	recneuron = 0;
	recstart = 53000 * 1000;
	recstop = 54000 * 1000;

	timestart = clock();

	// Model Loop
	for(step=1; step<=modsteps; step++) {
		//ttime = ttime + hstep;
		//neurotime = neurotime + hstep;
		ttime++;
		neurotime++;
		//if(countflag && (int)neurotime % (runtime / 100) == 0) netmod->mod->dispbox->SetCount(floor(neurotime)/modsteps*100);

		if(countflag && (int)neurotime % runtime100 == 0) {
			progevent.SetInt(floor(neurotime)/modsteps*100);  
			netbox->GetEventHandler()->AddPendingEvent(progevent);
		}
		
		if(step % 1000 == 0 && neuron->spikecount > 0) {
			plotevent.SetInt(floor(neurotime)/modsteps*100);  
			netbox->GetEventHandler()->AddPendingEvent(plotevent);
			if((*netmod->netflags)["realtime"]) Sleep(disprate);
		}

		// Osmo Net Sync
		if(osmomode && step % osmorate == 0) {
			while(step >= netmod->osmotime) {
				//netmod->diagmute->Lock();
				//diagbox->Write(text.Format("Neuron %d waiting step %d osmotime %d\n", neurodex, step, netmod->osmotime));
				//netmod->diagmute->Unlock();
				Sleep(100);
			}
			OsmoPress = netmod->OsmoStore[step / netmod->osmo_hstep];
			IrOsmoPress = (26 * (OsmoPress - 303)) / 1000;
			if(step % 100000 == 0) {
				//oxynetmod->diagmute->Lock();
				//oxynetmod->mod->diagbox->Write(text.Format("Neuron %d OsmoPress %.2f\n", neurodex, OsmoPress));
				//oxynetmod->diagmute->Unlock();
			}
		}


		// Signal Input     
		if(noiamp) noisig = noisig + (noimean - noisig) / noitau + noiamp * sqrt(hstep) * gaussian(0, 1); 
		if(signalmode) {
			synsig = noisig;
			epsprate1 = synsig / 1000; 
			ipsprate1 = epsprate1 * sigIratio; 	
		}
		else {
			synsig = psprate;
			epsprate1 = 0;
			ipsprate1 = 0;
		}

		
		// rate test
		//epsprate1 = 0;
		//ipsprate1 = 0;


		// PSP input signal
		nepsp = 0;
		nipsp = 0;
		nepsp1 = 0;
		nipsp1 = 0;
		nepsp2 = 0;

		if(inputgen) {
			nepsp = (neuron->dendinputE)[step];
			nipsp = (neuron->dendinputI)[step];
		}
		else {
			if(prototype == ramp) {
				//mod->diagbox->Write("set ramp\n"); 
				if(step < rampstart) rampinput = rampbase;
				if(step >= rampstart && step < rampstop) rampinput = rampinit  + (step - rampstart) * rampstep;  // * hstep
				if(step >= rampstop) rampinput = rampafter;
				if(rampinput < 0) rampinput = 0;
				epsprate = rampinput / 1000;
				synsig = rampinput;
			}
			if(prototype == rampcurve) {
				if(step < rampstart) rampinput = rampbase;
				if(step >= rampstart && step < rampstop) rampinput = rampinit  + rampmax - rampmax * exp(-rampgrad * (step-rampstart));
				if(step >= rampstop) rampinput = rampafter;
				if(rampinput < 0) rampinput = 0;
				epsprate = rampinput / 1000;
				synsig = rampinput;
			}

			// psp rate can be affected by randomly arriving psp but also by stimuli
			//totalepsprate = (epsprate + IrOsmoPress) * synvar;
			//totalipsprate = (ipsprate + IrOsmoPress) * pspRatio * synvar;

			totalepsprate = epsprate * synvar;
			totalipsprate = epsprate * pspRatio * synvar;

			if(totalepsprate > 0) {
				while(epspt < hstep) {
					erand = unif01(randgen);
					//erand = para_mrand01(neurodex);
					//erand = sfmt_genrand_real2(&sfmt);
					nepsp++;
					//epspt = -log(1 - para_mrand01(neurodex)) / totalepsprate + epspt;
					epspt = -log(1 - erand) / totalepsprate + epspt;
					//epspt = -log(1 - dis(randmt)) / totalepsprate + epspt;
				}
				epspt = epspt - hstep;
			}

			if(!flagError && epspt > 1000) {
				mod->diagbox->Write(text.Format("epspt %.10f  erand %.10f\n", epspt, erand));
				flagError = true;
			}

			if(totalipsprate > 0) {
				while(ipspt < hstep) {
					irand = unif01(randgen);
					//irand = para_mrand01(neurodex);
					//irand = sfmt_genrand_real2(&sfmt);
					nipsp++;
					//ipspt = -log(1 - para_mrand01(neurodex)) / totalipsprate + ipspt;
					ipspt = -log(1 - irand) / totalipsprate + ipspt;
					//ipspt = -log(1 - dis(randmt)) / totalipsprate + ipspt;
				}
				ipspt = ipspt - hstep;
			}

			
			if(epsprate1 > 0) {
				while(epspt1 < hstep) {
					erand = unif01(randgen);
					//erand = para_mrand01(neurodex);
					nepsp1++;
					epspt1 = -log(1 - erand) / epsprate1 + epspt1;
				}
				epspt1 = epspt1 - hstep;
			}

			if(ipsprate1 > 0) {
				while(ipspt1 < hstep) {
					erand = unif01(randgen);
					//irand = para_mrand01(neurodex);
					nipsp1++;
					ipspt1 = -log(1 - irand) / ipsprate1 + ipspt1;
				}
				ipspt1 = ipspt1 - hstep;
			}

			if(epsprate2 > 0) {
				while(epspt2 < hstep) {
					erand = unif01(randgen);
					//erand = para_mrand01(neurodex);
					//erand = sfmt_genrand_real2(&sfmt);
					nepsp2++;
					epspt2 = -log(1 - erand) / epsprate2 + epspt2;
				}
				epspt2 = epspt2 - hstep;
			}
			
		}

		inputPSP = nepsp * epspmag - nipsp * ipspmag;
		inputPSP1 = nepsp1 * epspmag - nipsp1 * ipspmag;

		// Input dynamics
		if(epspmag2) {
			if(epspsynchflag) nepsp2 = nepsp;   // synchronous AMPA and NMDA EPSPs
			inputPSP2 = inputPSP2 - (inputPSP2 * tauPSP2) * hstep + nepsp2 * epspmag2;
		}

		// Spiking model
		
		//pspsig = pspsig + (inputPSP2 * tauPSP2 - pspsig * tauMem) * hstep + inputPSP + inputPSP1;
		pspsig = pspsig + (inputPSP2 * tauPSP2 - pspsig * tauMem) + inputPSP + inputPSP1;

		//tHAP = tHAP - (tHAP * tauHAP) * hstep;
		//tDAP = tDAP - (tDAP * tauDAP) * hstep;
		//tAHP = tAHP - (tAHP * tauAHP) * hstep;
		//tAHP2 = tAHP2 - (tAHP2 * tauAHP2) * hstep;     // currently redundant hstep removed for speed optimization

		tHAP = tHAP - (tHAP * tauHAP);
		tDAP = tDAP - (tDAP * tauDAP);
		tAHP = tAHP - (tAHP * tauAHP);
		tAHP2 = tAHP2 - (tAHP2 * tauAHP2);

		tCa = tCa - (tCa - Ca_rest) * tauCa;
		tDyno = tDyno - tDyno * tauDyno;

		//tdendCa = tdendCa - hstep * tdendCa * taudendCa;
		tdendCa = tdendCa - tdendCa * taudendCa;
		storeDyno = storeDyno + kstoreDyno * tdendCa;  // - hstep * neuron->storeDyno / taustoreDyno;
		if(storeDyno > 10) storeDyno = 10;   

		// Osmosensitive Depolarisation
		//inputOsmo = Osmo * gOsmo;
		inputOsmo = gOsmo;

		// IKleak

		KLact = vox_tanh((tCa - Ca_rest - tDyno) / ka);  
		//IKL = gKL * (1 - KLact);
		IKL = gKL - gKL * KLact;

		V = Vrest + pspsig + inputOsmo - tHAP - tAHP - tAHP2 + tDAP - IKL; 


		// Secretion model

		if(netmod->secmode) {
			// Secretion dynamics
			//tB = tB - (tB * tauB) * hstep;   // broadening
			//tE = tE - (tE * tauE) * hstep;   // fast Ca2+
			//tC = tC - (tC * tauC) * hstep;   // slow Ca2+
			tB = tB - (tB * tauB);   // broadening
			tE = tE - (tE * tauE);   // fast Ca2+
			tC = tC - (tC * tauC);   // slow Ca2+

			//tC += - shstep * tC * tauC;          // slower Ca accumulation
			//tE += - shstep * tE * tauE;          // fast Ca for exocytosis

			// Inhibitory feedback from Ca currents in the axonal projection
			// For oxytocin there is a smaller negative feedback for the fast and slow Ca++
			//Cinh = 1 - pow(tC, Cgradient) / (pow(tC, Cgradient) + pow(Cthresh, Cgradient));		
			//Einh = 1 - pow(tE, Egradient) / (pow(tE, Egradient) + pow(Ethresh, Egradient)); 
			
			EKpow = tE * tE * tE * tE * tE;
			//Einh = 1 - EKpow / (EKpow + Ethresh * Ethresh * Ethresh * Ethresh * Ethresh);
			Einh = 1 - EKpow / (EKpow + Ethpow);
			CKpow = tC * tC * tC;  // * tC * tC;
			//Cinh = 1 - CKpow / (CKpow + Cthresh * Cthresh * Cthresh); // * Cthresh * Cthresh);
			Cinh = 1 - CKpow / (CKpow + Cthpow); 

			CaEnt = Einh * Cinh * (tB + Bbase); // Ca Entry		
			
			//secX = pow(tE, secExp) * alpha * tP;			// Rate of secretion (vesicle exocytosis)
			if(secExp == 3) secX = tE * tE * tE * alpha * tP;             // fixed vaso
			if(secExp == 2) secX = tE * tE * alpha * tP;                    // fixed oxy

			// Reserve Store (tR) and Releasable Pool (tP)
			if(tP < Pmax) fillP = beta * tR / Rmax; 
			else fillP = 0;

			tP = tP - secX + fillP;

			secBinX += secX;
		}

		// Plasma model

		if(netmod->secmode && netmod->plasmamode) {

			secXbuffer[buffdex++] = secX / netmod->numneurons;

			/*if(step >= 2000000 && step < 2000010) {
			oxynetmod->diagmute->Lock();
			oxynetmod->mod->diagbox->Write(text.Format("NeuroMod secX %.2f step %d secXtime %d\n", secX, step, magpop->secXtime));
			oxynetmod->diagmute->Unlock();
			}*/

			// New Buffered Secretion Rate code

			if(buffrate && step >= buffrate && step % buffrate == 0) {
				netmod->secmute->Lock();
				//for(i=0; i<buffrate; i++) magpop->secX[step - buffrate + i] += secXbuffer[i];
				for(i=0; i<buffrate; i++) secXpop[(step - buffrate + i) / plasma_hstep] += secXbuffer[i];
				magpop->secXcount[step / buffrate]++;
				if(magpop->secXcount[step / buffrate] == netmod->numneurons) magpop->secXtime = step;
				netmod->secmute->Unlock();
				buffdex = 0;
				if(neurodex == 0 && step < 2000) {
					//netmod->diagmute->Lock();
					netmod->mod->DiagWrite(text.Format("Neuron 0 buffer fill step %d secXtime %d secXpop %d buffer %d\n", step, magpop->secXtime, (step - buffrate)/plasma_hstep, buffrate));
					//netmod->diagmute->Unlock();
				}
			}

			// old non-buffered code removed here - see archived versions 2/5/22

			// bin recording of secretion rate and plasma concentration
			secRate1s += secX;
			plasmaRate1s += tOxyPlasma;		
			secRate60s += secX;
			secRate600s += secX;

			if((step % 1000) == 0) {
				neuron->Secretion[step/1000] = secRate1s; 
				//neuron->OxyPlasma[step/1000] = plasmaRate1s;
				secRate1s = 0;
				plasmaRate1s = 0;
			}

			if((step % 60000) == 0) {
				neuron->secLong[step/60000] = secRate60s * 60 / 1000;  // convert pg/min to ng/h
				secRate60s = 0;
			}

			// currently set for 10 minute window
			if((step % 600000) == 0) {
				neuron->secHour[step/600000] = secRate600s * 6 / 1000;  // convert pg/min to ng/h
				secRate600s = 0;
			}
		}  

		// Synthesis model

		// 'shstep' is 1 second scale time step

		// Synthesis V8               14th July 2021

		stimTS += (kTS * 0.001 * (tCa - Ca_rest) - stimTS * tauTS) * shstep;

		stimTL += (kTL * 0.001 * (tCa - Ca_rest) - stimTL * tauTL) * shstep;

		synthrate = (stimTL + basalTL) * mRNAstore;

		if(decaymode) mRNAstore += (synscale * stimTS - mRNAstore * mRNAtau) * shstep;	
		else mRNAstore += synscale * (stimTS - synthrate) * shstep;	
		
		if(mRNAmax && mRNAstore > mRNAmax) mRNAstore = mRNAmax;  

		//fillR = rateSR * synthrate; 
		//fillR = rateSR * synthrate * 0.00003;
		if(!synthdel) fillR = rateSR * synthrate * 0.001 * 0.03;
		else {
			if(step/synthrecrate >= synthdel) fillR = rateSR * synthrec[step/synthrecrate - synthdel] * 0.001 * 0.03;
			else fillR = rateSR * synthrec[0] * 0.001 * 0.03;
		}

		tR = tR + fillR - fillP;    // Reserve store - linking synthesis to secretion model

		if(step%1000 == 0 && step/1000 < magpop->maxtime) neuron->store[step/datsample] = tR;
		
		// Neuron Monitor
		if(neurodex == 0 && step%100 == 0 && step<1000000) {
			magpop->inputsignal[step/100] = synsig;	         // input signal recording
		}

		if(neurodex == 0 && step < 1000) {
			//fprintf(tofp, "Step %d	Input %.2f  EPSP rate %.2f  IPSP rate %.2f  InputSyn %.2f  HAP %.2f  V %.2f\n", 
			//	step, inputPSP, totalepsprate, totalipsprate, pspsig, HAP, V); 
		}

		if(neurodex == 0 && step < 1000000) {
			netmod->neurodata->pspsig[step] = pspsig;
			//netmod->oxyneurodata->synsig[step] = 10;
		}
		

		// Record - minute timescale, sampled
		if(step < 60000 * maxtimeLong) {
			if(step % 60000 == 0) { 
				//neuron->vasoLong[step/60000] = nsV;
				neuron->storeLong[step/60000] = tR;
				neuron->transLong[step/60000] = stimTS;
				//neuron->CaLong[step/60000] = tCa;
				//diagbox->Write(text.Format("step %d tR %.2f\n", step, tR));
				neuron->synthstoreLong[step/60000] = mRNAstore;
				neuron->synthrateLong[step/60000] = fillR * 3600;   // convert from pg/ms to ng/h
				if(neurodex == 0) magpop->inputLong[step/60000] = synsig;
			}
		}

		if(monitor && neurodex == 0) {
			if(step % datsample == 0 && step < 1000000 * datsample) {
				neurorecord->stimTS[step/datsample] = stimTS;
				neurorecord->stimTL[step/datsample] = stimTL;
				neurorecord->mRNAstore[step/datsample] = mRNAstore;
				neurorecord->Ca[step/datsample] = tCa;
			}
		}

		// Synth delay store
		//
		// record synthesis rate every minute for delay recall

		if(step%synthrecrate == 0) { 
			synthdex = step/synthrecrate;
			if(synthdex > 100000) synthdex = synthdex % 100000;
			synthrec[synthdex] = synthrate;
			//diagbox->Write(text.Format("neuromod synthrec step %d dex %d data %.4f\n", step, synthdex, synthrate));
		}


		/*if(monitor && neurodex == recneuron)            // 15/7/20, diagnostics used to find random number error
			if(step >= recstart && step < recstop) {
				neurorecord->V[step - recstart] = V;
				neurorecord->syn[step - recstart] = pspsig;
				neurorecord->psp[step - recstart] = epspt; // erand; //epspt; //nipsp;   // inputPSP;
				neurorecord->rand[step - recstart] = erand; 
			}*/

		// Spiking
		if(V > Vthresh && ttime >= absref) {

			// record spike time
			if(neuron->spikecount < maxspikes) {
				neuron->times[neuron->spikecount] = neurotime;
				neuron->spikecount++;
			}

			neuron->spikecount2++;

			// Spike incremented variables

			// Calcium
			tCa = tCa + kCa;

			// afterpotentials
			tAHP = tAHP + kAHP;
			tHAP = tHAP + kHAP;
			tDAP = tDAP + kDAP;

			if(AHP2mode) {
				if(tCa >= aAHP2) tAHP2 = tAHP2 + kAHP2 * (tCa - aAHP2);
			}
			else tAHP2 = tAHP2 + kAHP2;

			// secretion variables
			if(netmod->secmode) {
				tB = tB + kB;
				tE = tE + kE * CaEnt;
				tC = tC + kC * CaEnt;			
			}

			// Dynorphin
			if(dynostoreflag) {
				if(storeDyno > spikeDyno) {
					tDyno = tDyno + kDyno;
					storeDyno = storeDyno - spikeDyno;
				}
			}
			else tDyno = tDyno + kDyno;
		}
	}


	// Store final mRNA store and reserve store value for sequential runs
	if(!neuron->netinit) {
		neuron->mRNAinit = mRNAstore;
		(*neuron->synthparams)["mRNAinit"] = mRNAstore;
	}

	if(!neuron->storereset) {
		(*neuron->secparams)["Rinit"] = tR;
		neuron->storeinit = tR;
	}

	//oxynetmod->diagmute->Lock();
	//oxynetmod->mod->diagbox->Write(text.Format("Neuron %d finished step %d secXtime %d\n", neurodex, step, magpop->secXtime));
	//oxynetmod->diagmute->Unlock();
	
	//fclose(tofp);
	delete [] secXbuffer;
	delete [] synthrec;
}



float calcLognorm(float mean, float StDev)
{
	float lognMean, lognSD;
	float lognresult; 

	// lognormal distributed numbers
	std::mt19937 randgen;
	randgen.seed(time(NULL));

	//// Mean and SD need to be recalculated. 
	//// https: //stackoverflow.com/questions/23699738/how-to-create-a-random-number-following-a-lognormal-distribution-in-excel
	lognMean = log(pow(mean, 2)/sqrt(pow(mean, 2) + pow(StDev, 2)));
	lognSD = sqrt (log( (pow(mean, 2) + pow(StDev, 2)) / pow(mean, 2) ));

	////the parameters are (mean, standard deviation) - default is (0.0, 1.0)
	std::lognormal_distribution<float> lognorm(lognMean, lognSD);

	lognresult = lognorm(randgen);

	return lognresult;
}
