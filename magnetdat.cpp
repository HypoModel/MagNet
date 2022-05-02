/*
*  magnetdat.cpp
*  HypoModel
*
*  Created by Duncan MacGregor
*  June 2017
*
*/


#include "magnetmodel.h"


MagNeuron::MagNeuron()
{
	spikeparams = new ParamStore();
	secparams = new ParamStore();
	sigparams = new ParamStore();
	dendparams = new ParamStore();
	synthparams = new ParamStore();
	protoparams = new ParamStore();
	spikecount = 0;
	
	maxtimeRate1s = 100000;
	maxtimeRate10s = 10000;
	maxtime = maxtimeRate1s;
	maxtimeLong = 35000;

	Secretion.setsize(maxtimeRate1s);

	store.setsize(maxtime); //, mainwin->diagbox->textbox, "b");
	storeLong.setsize(maxtimeLong); //, mainwin->diagbox->textbox, "storeLong"); 
	transLong.setsize(maxtimeLong); // mainwin->diagbox->textbox, "transLong"); 
	//CaLong.setsize(maxtimeLong); //, mainwin->diagbox->textbox, "CaLong"); 
	synthstoreLong.setsize(maxtimeLong); //, mainwin->diagbox->textbox, "synthstoreLong"); 
	synthrateLong.setsize(maxtimeLong); 
	secLong.setsize(maxtimeLong);
	secHour.setsize(maxtimeLong);

	initflag = false;
	synvar = 1;
	mRNAinit = 0;
	storeinit = 0;
}


MagNeuron::~MagNeuron()
{
	delete spikeparams;
	delete secparams;
	delete sigparams;
	delete dendparams;
	delete synthparams;
	delete protoparams;
}


void MagNeuron::StoreClear()
{
	Secretion.reset();
	store.reset();
	storeLong.reset();
	CaLong.reset();
	synthstoreLong.reset();
	transLong.reset();
}


MagPop::MagPop()
{
	maxtime = 200000;
	maxtimeRate1s = maxtime;
	maxtimeRate10ms = maxtime * 100;
	maxtimeRate1ms = maxtime * 1000;
	maxtimeLong = 35000;  // minutes, 30000 sufficient for 20 days simulation and recording

	OxySecretionNet.setsize(maxtimeRate1s);
	OxyPlasmaNet.setsize(maxtimeRate1s);
	PlasmaNaConc.setsize(maxtimeRate1s);
	EVFNaConc.setsize(maxtimeRate1s);
	DiffNaGrad.setsize(maxtimeRate1s);
	ICFGrad.setsize(maxtimeRate1s);
	ICFVol.setsize(maxtimeRate1s);
	EVFNaVol.setsize(maxtimeRate1s);

	NetSecretion4s.setsize(maxtimeRate1s);

	evfNaConcTemp.setsize(maxtimeRate1ms);   // evfNaConc variable buffer for calculating moving 2s average

	EVFNaConc.setsize(maxtimeRate10ms);
	PlasmaNaConc.setsize(maxtimeRate10ms);
	OsmoPress1s.setsize(maxtimeRate10ms);
	inputsignal.setsize(maxtimeRate10ms);
	netsignal.setsize(maxtimeRate10ms);

	inputLong.setsize(maxtimeLong);
	plasmaLong.setsize(maxtimeLong);
	netsecLong.setsize(maxtimeLong);
	netsecHour.setsize(maxtimeLong);
	secLong.setsize(maxtimeLong);
	secHour.setsize(maxtimeLong);
	transLong.setsize(maxtimeLong);

	secX.setsize(maxtimeRate1ms);
	secXcount.setsize(maxtime * 100);

	storeLong.setsize(maxtimeLong);
	synthstoreLong.setsize(maxtimeLong);
	storesum.setsize(maxtime);
	storesumLong.setsize(maxtimeLong);
	storesumNorm.setsize(maxtimeLong);
	synthstoresumLong.setsize(maxtimeLong);
	synthratesumLong.setsize(maxtimeLong);

	srate1s.setsize(maxtime);
	srate10s.setsize(maxtime/10);
	srate30s.setsize(maxtime/30);
	srate300s.setsize(maxtime/300);
	srate600s.setsize(maxtime/600);
}


void MagPop::StoreClear()
{
	OxySecretionNet.reset();
	OxyPlasmaNet.reset();
}


void MagPop::PopSum()
{
	int i, step, min;

	//std::vector<int> neuron_srate;
	//neuron_srate.resize(maxtime);

	// Clear store and rate counts
	//for(step=0; step<runtime && step<maxtime; step++) {
	for(step=0; step<maxtime; step++) {
		//vasosum[step] = 0;
		//srate[step] = 0;
		storesum[step] = 0;
	}

	for(min=0; min<maxtimeLong; min++) {
		storesumLong[min] = 0;
		synthstoresumLong[min] = 0;
		synthratesumLong[min] = 0;
	}
	numspikes = 0;

	storesum.reset();
	storesumLong.reset();
	storesumNorm.reset();
	synthstoresumLong.reset();
	synthratesumLong.reset();

	popfreq = 0;

	// Sum release and spike rate over population
	for(i=0; i<numneurons; i++) {
		//vasomod->currvaso->neurocalc(&(neurons[i]));
		//(*neurons)[i].ratecalc(neuron_srate);
		popfreq += (*neurons)[i].freq;
		for(step=0; step<runtime && step<maxtime; step++) {
			//vasosum[step] += neurons[i].vaso[step] / numcells;
			//vasomod->vasopop->srate[step] += vasomod->currvaso->srate[step];
			//srate[step] += (*neurons)[i].srate[step];
			//srate1s[step] += neuron_srate[step];
			storesum[step] += (*neurons)[i].store[step];
		}

		for(min=0; min<runtime/60; min++) {
			//vasosumLong[min] += neurons[i].vasoLong[min] / numcells;
			storesumLong[min] += (*neurons)[i].storeLong[min] / numneurons;
			storesumNorm[min] = storesumLong[min] / 10;
			synthstoresumLong[min] += (*neurons)[i].synthstoreLong[min] / numneurons;
			synthratesumLong[min] += (*neurons)[i].synthrateLong[min] / numneurons;
		}
	}

	popfreq = popfreq / numneurons;

	// Total spike sum
	//for(step=0; step<runtime; step++)
	//	numspikes += srate[step];
}


MagNeuroDat::MagNeuroDat()
{
	int storesize = 1000000;
	int sizeLong = 20000;
	int sizeSmall = 1000;

	//diagbox = main->diagbox;
	secP.setsize(sizeSmall, true);
	secR.setsize(sizeSmall, true);
	secX.setsize(sizeSmall, true);

	pspsig.setsize(storesize);
	//inputrate.setsize(storesize);
	Ca.setsize(storesize);
	V.setsize(storesize);
	syn.setsize(storesize);
	psp.setsize(storesize);
	rand.setsize(storesize);

	stimTL.setsize(storesize);
	stimTS.setsize(storesize);
	mRNAstore.setsize(storesize);
}


