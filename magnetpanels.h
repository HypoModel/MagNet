#pragma once


#include <hypomodel.h>
#include <magnetdat.h>


class MagNetModel;



class OxySignalBox: public ParamBox
{
public:
	MagNetModel *mod;

	OxySignalBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
	void OnRun(wxCommandEvent& event);
};


class MagNetProtoBox : public ParamBox
{
public:
	MagNetModel *mod;
	wxStaticText *currentinput;
	wxStaticText *currentpulse;
	wxStaticText *currentrange;
	wxStaticText *status;

	MagNetProtoBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);

	void OnRun(wxCommandEvent& event);
};


class MagNetBox: public ParamBox
{
public:
	MagNetModel *mod;
	wxCheckBox *seedcheck;
	wxCheckBox *initcheck;
	wxCheckBox *storecheck;

	MagNetBox(MagNetModel *mod, MainFrame *main, const wxString& title, const wxPoint& pos, const wxSize& size);
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
	MagNetModel *mod;
	wxCheckBox *synccheck;

	void OxyPanel();
	void VasoPanel();

	MagSpikeBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class MagSynthBox: public ParamBox
{
public:	
	MagSynthBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class OxyDendBox: public ParamBox
{
public:
	MagNetModel *mod;
	wxCheckBox *synccheck;

	OxyDendBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class OxySecBox: public ParamBox
{
public:
	MagNetModel *mod;
	wxCheckBox *synccheck;

	OxySecBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class OxyNeuroDataBox: public ParamBox
{
public:
	MagNetModel *mod;
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
	OxyNeuroDataBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
};


class MagNetGridBox : public GridBox
{
public:
	MagNetModel *mod;

	void OnPlot(wxCommandEvent& event);

	MagNetGridBox(MagNetModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size, int rows, int cols);
};