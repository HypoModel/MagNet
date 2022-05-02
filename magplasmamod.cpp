
/*
*  magplasmamod.cpp
*  HypoModel
*
*  Created by Duncan MacGregor.
*
*	
*
*
*/


#include "magnetmodel.h"


MagPlasmaMod::MagPlasmaMod(MagNetMod *oxynetmod)
	: wxThread(wxTHREAD_JOINABLE)
{
	netmod = oxynetmod;
	mod = netmod->mod;
	magpop = mod->magpop;

	ParamStore *secparams = mod->secbox->GetParams();
	ParamStore *secflags = mod->secbox->modflags;
	
	// Diffusion and Clearance Parameters
	halflifeOxyClear = (*secparams)["ClearHL"];
	halflifeOxyDiff = (*secparams)["DiffHL"];
	PlasmaVol = (*secparams)["VolPlasma"];
	EVFVol = (*secparams)["VolEVF"];
	plasma_hstep = (*secparams)["plasma_hstep"];

	diff_flag = (*secflags)["diff_flag"];
}


void *MagPlasmaMod::Entry()
{
	plasmamodel();
	return NULL;
}


void MagPlasmaMod::plasmamodel()
{
	int i, step;
	int runtime, modtime;
	wxString text;
	double plasmatime = 0;

	double DiffRate;
	double tauOxyClear, tauOxyDiff;
	double netsecRate1s, netplasmaRate1s;
	double netsecRate4s;
	double plasmaRate60s, netsecRate60s;
	double netsecRate1h;

	tauOxyClear = log((double)2) / (halflifeOxyClear * 1000);
	tauOxyDiff = log((double)2) / (halflifeOxyDiff * 1000);
	
	runtime = netmod->runtime * 1000;
	modsteps = runtime / plasma_hstep;
	buffrate = netmod->buffrate / plasma_hstep;

	if(!buffrate) return;

	netsecRate1s = 0;
	netsecRate4s = 0;
	netplasmaRate1s = 0;
	plasmaRate60s = 0;
	netsecRate60s = 0;
	netsecRate1h = 0;
	magpop->OxySecretionNet.reset();
	magpop->OxyPlasmaNet.reset();

    /*
	netmod->diagmute->Lock();
	netmod->mod->diagbox->Write(text.Format("PlasmaMod running secXtime %d modsteps %d\n", magpop->secXtime, modsteps));
	netmod->diagmute->Unlock();
     */
    mod->DiagWrite(text.Format("PlasmaMod running secXtime %d modsteps %d\n", magpop->secXtime, modsteps));
    

	// Model Loop
	for(step=1; step<=modsteps; step++) {
		plasmatime += plasma_hstep;
		// Wait for secX summation buffer
		if(step % buffrate == 0) {
			while(plasmatime > magpop->secXtime) {
				//oxynetmod->diagmute->Lock();
				//oxynetmod->mod->diagbox->Write(text.Format("PlasmaMod waiting step %d secXtime %d\n", step, oxypop->secXtime));
				//oxynetmod->diagmute->Unlock();
				Sleep(100);
			}
			/*if(step < 10000) {
				oxynetmod->diagmute->Lock();
				oxynetmod->mod->diagbox->Write(text.Format("PlasmaMod buffer scaling step %d secXtime %d\n", step, oxypop->secXtime));
				oxynetmod->diagmute->Unlock();
			}*/
		}
		
		// Diffusion Rate: will be positive or negative in one or another way depending on the oxytocin concentration in each compartment.
		if(!diff_flag) DiffRate = 0;
		else DiffRate = (netmod->tOxyPlasma / PlasmaVol - netmod->tOxyEVF / EVFVol) * (PlasmaVol + EVFVol) / 2; // the pressure is total amount, not from the amount/ml	

		// If [OxPlasma] > [OxEVF] -> {DiffRate > 0} -> tOxyPlasma will give plasma to tOxyEVF
		// If [OxPlasma] < [OxEVF] -> {DiffRate < 0} -> tOxyPlasma will receive plasma from tOxyEVF

		/*if(step >= 2000000 && step < 2000010) {
			oxynetmod->diagmute->Lock();
			oxynetmod->mod->diagbox->Write(text.Format("PlasmaMod secX %.4f step %d secXtime %d\n", oxypop->secX[step], step, oxypop->secXtime));
			oxynetmod->diagmute->Unlock();
		}*/

		netmod->tOxyPlasma = netmod->tOxyPlasma + plasma_hstep * (magpop->secX[step] - (netmod->tOxyPlasma * tauOxyClear + DiffRate * tauOxyDiff));  // Oxytocin Plasma Concentration
		netmod->tOxyEVF = netmod->tOxyEVF + plasma_hstep * (DiffRate * tauOxyDiff);

		//netsecRate1s =+ netmod->netsecX;
		netsecRate1s += magpop->secX[step];
		netsecRate4s += magpop->secX[step];
		netplasmaRate1s += netmod->tOxyPlasma;	
		plasmaRate60s += netmod->tOxyPlasma;
		netsecRate60s += magpop->secX[step];         
		netsecRate1h += magpop->secX[step];                  // long timescale secretion rate for fitting to Robinson 1989 

		if(step % (1000 / plasma_hstep) == 0) {
			magpop->OxySecretionNet[step/(1000/plasma_hstep)] = mod->popscale * netsecRate1s; 
			//netmod->mod->oxypop->OxySecretionNet[step/1000] = 10;
			magpop->OxyPlasmaNet[step/(1000/plasma_hstep)] = (netplasmaRate1s / 1000) / PlasmaVol;
			netsecRate1s = 0;
			netplasmaRate1s = 0;
		}

		if(step % (4000 / plasma_hstep) == 0) {
			magpop->NetSecretion4s[step/(4000/plasma_hstep)] = mod->popscale * netsecRate4s;
			netsecRate4s = 0;
		}

		if(step % (60000 / plasma_hstep) == 0) {
			magpop->plasmaLong[step/(60000/plasma_hstep)] = (plasmaRate60s / 60000) / PlasmaVol;
			magpop->netsecLong[step/(60000/plasma_hstep)] = mod->popscale * netsecRate60s * 60 / 1000;    // convert pg/min to ng/h
			plasmaRate60s = 0;
			netsecRate60s = 0;
		}

		// notional hour rate, currently set to 10 minute
		if(step % (600000 / plasma_hstep) == 0) {
			magpop->netsecHour[step/(600000/plasma_hstep)] = mod->popscale * netsecRate1h * 6 / 1000; 
			netsecRate1h = 0;
		}
	}

    /*
	netmod->diagmute->Lock();
	netmod->mod->diagbox->Write(text.Format("PlasmaMod finished secXtime %d plasma maxdex %d\n", magpop->secXtime, magpop->OxyPlasmaNet.maxdex()));
	netmod->diagmute->Unlock();
     */
    
    mod->DiagWrite(text.Format("PlasmaMod finished secXtime %d plasma maxdex %d\n", magpop->secXtime, magpop->OxyPlasmaNet.maxdex()));
}
