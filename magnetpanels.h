/*
*  magnetpanels.h
*  
*  Created by Duncan MacGregor
*  University of Edinburgh 2022
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*
*/



#pragma once


#include "hypomain.h"
#include "magnetdat.h"


class MagNetMod;


class MagGenBox: public ParamBox
{
public:	
	wxString gentags[20];
	int numgen;
	MagNetMod *mod;

	MagGenBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
	void OnGenerate(wxCommandEvent& event);
	void OnZero(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
};


class MagSignalBox: public ParamBox
{
public:
	MagNetMod *mod;

	MagSignalBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
	void OnRun(wxCommandEvent& event);
};


class MagNetProtoBox : public ParamBox
{
public:
	MagNetMod *mod;
	wxStaticText *currentinput;
	wxStaticText *currentpulse;
	wxStaticText *currentrange;
	wxStaticText *status;

	MagNetProtoBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size);

	void OnRun(wxCommandEvent& event);
};


class MagNetBox: public ParamBox
{
public:
	MagNetMod *mod;
	wxCheckBox *seedcheck;
	wxCheckBox *initcheck;
	wxCheckBox *storecheck;

	MagNetBox(MagNetMod *mod, MainFrame *main, const wxString& title, const wxPoint& pos, const wxSize& size);
	void OnParamStore(wxCommandEvent& event);
	void OnParamLoad(wxCommandEvent& event);
	//void OnBox(wxCommandEvent& event);

	void OnRun(wxCommandEvent& event);
	void OnProgress(wxCommandEvent& event);
	void OnPlot(wxCommandEvent& event);
	void NetStore();
	void NetLoad();
};


class MagSpikeBox: public ParamBox
{
public:
	MagNetMod *mod;
	wxCheckBox *synccheck;

	void OxyPanel();
	void VasoPanel();
	//void NeuroGen(ParamStore *);

	MagSpikeBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class MagSynthBox: public ParamBox
{
public:	
	MagSynthBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class MagDendBox: public ParamBox
{
public:
	MagNetMod *mod;
	wxCheckBox *synccheck;

	MagDendBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class MagSecBox: public ParamBox
{
public:
	MagNetMod *mod;
	wxCheckBox *synccheck;

	MagSecBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class MagNeuroDataBox: public ParamBox
{
public:
	MagNetMod *mod;
	wxCheckBox *synccheck;
	int neurodex; // index for going through all the neurones of the network
	int neurocount;

	//MagNeuron *neurons;

	wxTextCtrl *datneuron;
	wxSpinButton *datspin;

	void NeuroData();
	void PanelData(NeuroDat *data);
	void OnNext(wxSpinEvent& event);
	void OnPrev(wxSpinEvent& event);
	void OnEnter(wxCommandEvent& event);
	MagNeuroDataBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class MagNetGridBox : public GridBox
{
public:
	MagNetMod *mod;

	void OnPlot(wxCommandEvent& event);

	MagNetGridBox(MagNetMod *mod, const wxString& title, const wxPoint& pos, const wxSize& size, int rows, int cols);
};
