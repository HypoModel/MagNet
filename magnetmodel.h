/*
*  magnetmodel.h
*
*  Created by Duncan MacGregor.  
*
*	Classes:
*		- "MagNeuron: public NeuroDat"    --->   Single neuron parameters and variables (spiking ones from class NeuroData, the rest included here)
*		- "MagNetDat"    --->   Just getting parameters for the Network (see magnetdat.cpp)
*		- Boxes for the network and the single neuron starting parameters  (see magnetpanels.cpp)
*		- "MagNeuroMod : public wxThread"   --->  Thread for working with a single neuron  (see magneuromod.cpp)
*		- "MagNetMod : public ModThread"   --->  class to work with threats, in this case with the single neuron threads defined by MagNeuroMod  (see magnetmod.cpp)
*		- "MagNetModel : public Model"   --->  to run the network threads and coordinate the graphs, boxes and data asociated  (see magnetmodel.cpp)
*		
*/


#ifndef MAGNETMOD_H
#define MAGNETMOD_H


#include "wx/wx.h"
#include <hypomodel.h>
#include <magnetdat.h>
#include <magnetpanels.h>
#include <hyponeuro.h>


enum {
	ID_oxynetmflag = 9000,
	ID_singletrans,
	ID_netgen,
	ID_seedgen,
	ID_setactive,
	ID_netdiag,
	ID_ipInfusion,
	ID_ivInfusion,
	ID_secmode,
	ID_plasmamode,
	ID_AHP2mode
};

class MagNetModel;
class MagNetMod;

class MagPlasmaMod : public wxThread
{
public:
	MagNetMod *netmod;
	MagNetModel *mod;
	MagPop *magpop;

	int modsteps;
	int netrate;
	int buffrate;
	int plasma_hstep;
	bool diff_flag;

	double PlasmaVol, EVFVol;
	double halflifeOxyClear, halflifeOxyDiff;

	MagPlasmaMod(MagNetMod *magnetmod);
	virtual void *Entry();

	void plasmamodel();
};


// Neuron model thread class
class MagNeuroMod : public wxThread
{
public:
	MagNeuron *neuron;  // neuron object
	MagNetMod *netmod;  // parent network model thread
	MagPop *magpop;     // population data store
	int neurodex;  // neuron index
	MagSpikeBox *spikebox;   // spiking model parameter box
	MagNetBox *netbox;
	MagNetModel *mod;
	DiagBox *diagbox;
	MagNeuroDat *neurorecord;

	int maxtime;
	int maxtimeLong;
	int prototype;

	// Model time steps
	int modsteps;
	double hstep, shstep, syn_hstep;

	// Spiking Parameters
	double Vrest, Vthresh;
	double pspmag, psprate, iratio;
	double kHAP, halflifeHAP;
	double kDAP, halflifeDAP;
	double halflifeMem;
	double kAHP, halflifeAHP;
	double synvar;

	// Vaso Spiking Parameters
	double kAHP2, halflifeAHP2, aAHP2;
	double kDyno, halflifeDyno;
	double kCa, halflifeCa, Ca_rest;
	double gKL, ka;
	double gOsmo;

	// Dynamic Dynorphin Parameters
	double kstoreDyno, halflifestoreDyno; 
	double spikeDyno; 
	double tauDynoup;
	double kdendCa, halflifedendCa; 

	// Secretion Parameters
	double Rmax, Rinit, Pmax;
	double kB, halflifeB, Bbase;
	double kE, halflifeE;
	double kC, halflifeC;
	double Cgradient, Cthresh;
	double Egradient, Ethresh;
	double secExp;
	double alpha, beta;

	// Synthesis Parameters         
	double stimTS, stimTL;
	double synthrate, mRNAstore;
	double fillP, fillR;

	double kTL, halflifeTL, tauTL;
	double basalTL, maxTL;
	double kTS, halflifeTS, tauTS;
	double mRNAinit, mRNAmax;
	double mRNAtau, mRNAhalflife;
	int synthdel; 
	double rateSR;  // scales synthesis units to reserve store units 
	double synscale;
	double vsynrate, vtrans;    // still in use 29/6/20?

	// Synthesis Flags
	int transflag;
	int transmode, synthmode;
	int scalemode;
	int polymode;
	int decaymode;

	// Protocol Parameters
	double rampbase;
	double rampstart, rampstop;
	double rampinit, rampstep;
	double rampafter;
	double rampmax, rampgrad;    // new 4/8/21 parameters for ramp curved to limit (max)


	double PlasmaVol, EVFVol;
	double halflifeClear, halflifeDiff;
	double BasalNaConc;
	int plasma_hstep;

	int netrate, osmorate;
	int buffrate;
	int osmomode;
	int disprate;
	unsigned long modseed;

	// NMDA synapse EPSPs
	double pspmag2, psprate2;
	double halflifePSP2;

	// Flags
	bool ivInfusionflag;
	bool ipInfusionflag;
	bool epspsynchflag;
	int modmode;    // switch between Oxy and Vaso parameter sets, oxy=0, vaso=1
	int AHP2mode;   // 0 for basic, 1 for Ca-threshold (PLoS 2012)
	bool dynostoreflag;  // dynamic store based dynorphin (PLoS 2012)
	
	// Noise Signal
	double noimean, noitau, noiamp;
	double sigIratio;
	bool noisemode, signalmode;

	MagNeuroMod(int index, MagNeuron *neuron, MagNetMod *oxynetmod);

	// running the model for a single neuron (each time)
	void neuromod();
	virtual void *Entry();
	//void calcLognorm();
};


// Class to manage all the threads of single neuron calculations. 
class MagNetMod : public ModThread
{
public:
	std::vector<MagNeuron> &neurons;
	MagNetModel *mod;
	MagNetDat *netdata;
	MagNeuroMod *neurothread[1000];  // for storaging data up to 1000 single neurones
	MagSpikeBox *spikebox;
	MagSynthBox *synthbox;
	MagNetBox *netbox;
	MagNeuroDat *neurodata;
	MagPlasmaMod *plasmathread;
	//OxyOsmoMod *osmothread;
	//OxySigMod *sigthread;
	MagPop *magpop;

	ParamStore *netflags;
	ParamStore *netparams;

	wxMutex *secmute;
	wxMutex *osmomute;
	int netrate, osmorate, buffrate;
	int osmo_hstep;
	int secmode, osmomode, plasmamode;
	unsigned long modseed;

	double netsecX;
	double netsecRate1s, netplasmaRate1s;
	double tPlasma, tEVF;
	double OsmoPress;  // not currently used, see OsmoStore

	datdouble OsmoStore;   // osmotic pressure buffer for feeding neuron threads
	int osmotime;

	// Protocol Flags
	bool rampflag;
	int prototype;

	//Protocol Parameters
	int *rampstart, *rampstop;
	double *rampbase, *rampinit;
	double *rampstep, *rampinput;
	double *rampafter;

	int runtime;
	int numneurons;
	wxMutex *diagmute;
	bool initflag;

	void Initialise();
	void RunNet();
	void Export2file(int, wxString, datdouble);
	//void NeuroGen();    // moved to MagNetModel
	int InputGen();
	void SecretionAnalysis();
	void RunRange();

	MagNetMod(MagNetModel *mod);
	virtual void *Entry();
};


class MagNetModel : public NeuroMod
{
public:
	MagNetBox *netbox;
	MagSpikeBox *spikebox;
	MagSecBox *secbox;
	MagNeuroDataBox *neurodatabox;
	MagSignalBox *signalbox;
	MagDendBox *dendbox;
	MagNetProtoBox *protobox;
	MagSynthBox *synthbox;

	MagNetDat *netdata;
	MagNeuroDat *neurodata;

	ParamBox *dispbox;

	std::vector<MagNeuron> modneurons;
	SpikeDat *currmodneuron;  // currneuron is a pointer to a single instance of the class. 
	SpikeDat *netdat;
	SpikeDat *netneuron;

	// Protocol Data Plots
	datdouble rangeref;
	datdouble rangedata[10];
	int rangeindex, rangecount, rangedatamax;

	int neuroindex;
	int celltypes;
	MagPop *magpop;

	datdouble nethist;
	bool netready;

	int numneurons;
	double popscale;
	int datsample;

	MagNetModel(int, wxString, HypoMain *);
	~MagNetModel();

	void RunModel();  
	void GraphData();
	void GSwitch(GraphDisp *gpos, ParamStore *gflags);
	void StoreClear();
	void DataOutput();
	void ScaleConsoleAbove(ScaleBox *, int condex);
	void ScaleConsoleBelow(ScaleBox *, int condex);
	//void DataSelect(double, double);
	//int SoundLink(SpikeDat **, datdouble **);
	int SoundLink(SoundBox *);
	void SoundOn();
	void RangePlot(TextGrid *textgrid);
	void ModClose();
    void OnModThreadCompletion(wxThreadEvent&);
    void OnEndRun(wxCommandEvent&);
	void ParamScan();
	void NeuroGen();   // moved from MagNetMod
};


#endif

