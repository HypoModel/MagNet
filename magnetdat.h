

#ifndef MAGNETDAT_H
#define MAGNETDAT_H

#include <hypomodel.h>


// 'MagNeuron' single magnocellular neuron class derived from NeuroDat
// NeuroDat contains spiking model variables and analysis storage for FR, ISI and hazard.
//
class MagNeuron : public NeuroDat
{
public:
	int active;
	int setactive;
	int index;
	int maxtime;
	int maxtimeLong;
	int maxtimeRate1s, maxtimeRate10s;

	ParamStore *spikeparams;    // Used for heterogeneity to store independent spiking model parameter set
	ParamStore *secparams;      // Not currently used
	ParamStore *sigparams;      // Input signal parameters
	ParamStore *dendparams; 
	ParamStore *synthparams;
	ParamStore *protoparams;

	// Pre-generated PSP counts for non-independent network inputs 
	unsigned char *dendinputE;
	unsigned char *dendinputI;

	// secretion and diffusion variable recording arrays
	datdouble Secretion;
	datdouble Plasma;
	datdouble secLong;    //  1 minute timescale
	datdouble secHour;    //  10 minute or 1 hour timescale

	// Stored initial values
	double synvar;
	double mRNAinit, storeinit;
	bool initflag, netinit;
	bool storereset;

	datdouble store;    // reserve store
	datdouble storeLong;
	datdouble transLong;
	datdouble CaLong;
	datdouble synthstoreLong;
	datdouble synthrateLong;

	MagNeuron();
	~MagNeuron();
	void StoreClear();
};


// 'MagPop' oxytocin or vasopressin population/network storage class
//
// 'neurons' contains MagNeuron array
// MagPop() does post run analysis and population data summation
// SpikeDat and datdouble pointers contain data for currently selected neuron
//
class MagPop{
public:
	int runtime;
	int numneurons;
	int numspikes;
	double oxymean;
	double oxytotal;
	std::vector<MagNeuron> *neurons;
	double popfreq;
	double popsd;
	double ratemean;
	double rateSD;
	SpikeDat *currneuron;

	int maxtime;
	int maxtimeRate1s;
	int maxtimeRate10ms;
	int maxtimeRate1ms;
	int maxtimeLong;

	datdouble *OxySecretion;     // Pointers to current neuron's secretion and plasma records
	datdouble *OxyPlasma;

	//datdouble OxyOsmo;
	datdouble OxySecretionNet;
	datdouble NetSecretion4s;
	datdouble OxyPlasmaNet;
	datdouble inputsignal;
	datdouble netsignal;
	datdouble inputLong;
	datdouble plasmaLong;
	datdouble netsecLong;
	datdouble netsecHour;
	datdouble secLong;
	datdouble secHour;
	datdouble transLong;

	// Synthesis and Stores
	datdouble storeLong;
	datdouble synthstoreLong;
	datdouble synthrateLong;
	datdouble storesum;
	datdouble storesumLong;
	datdouble storesumNorm;
	datdouble synthstoresumLong;
	datdouble synthratesumLong;

	// Osmotic Pressure model
	// 1 second bins
	datdouble PlasmaNaConc;
	datdouble EVFNaConc;
	datdouble DiffNaGrad;
	datdouble ICFGrad;
	datdouble ICFVol;
	datdouble EVFNaVol;
	datdouble OsmoPress1s;

	double PlasmaNaConcTemp;
	//double OsmoPress;
	int osmotime;
	//int osmostep;

	// Analysis
	double secmean, secmean_4s;
	double secIoD, secIoD_4s;

	// Summed spike rate
	datdouble srate1s, srate10s;
	datdouble srate30s;
	datdouble srate300s, srate600s;
	
	// 1ms bin
	datdouble evfNaConcTemp;

	// Summed Population secretion rate
	datdouble secX;
	datint secXcount;
	int secXtime;

	MagPop();
	//void Output(wxString tag);
	void PopSum();
	void StoreClear();
};


// 'MagNeuroDat' data storage class with recording model variables during a run
//
class MagNeuroDat{
public:
	int rectime;
	int datsample;
	
	SpikeDat *spikedat;
	Model *mod;
	DiagBox *diagbox;

	datdouble secX;
	datdouble secP;
	datdouble secR;

	datdouble Ca;
	datdouble V;
	datdouble syn;
	datdouble psp;
	datdouble rand;

	datdouble pspsig;
	//datdouble inputrate;

	datdouble stimTS;
	datdouble stimTL;
	datdouble mRNAstore;
	datdouble synthratestore;
	
	MagNeuroDat();
};

// 'MagNetDat' network/population class containing network parameters and neuron array link
//
class MagNetDat{
public:
	int runtime;
	int numneurons;

	MagNeuron *neurons;
	SpikeDat *currneuron;

	int maxtime;
	int maxtimeLong;

	MagNetDat();
};



#endif