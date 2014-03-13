/*
 * Fit.h
 *
 *  Created on: Dec 24, 2013
 *      Author: philip
 */

#ifndef FIT_H_
#define FIT_H_

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"
#include <TMinuit.h>
#include <TMath.h>
#include "TObject.h"
#include "TROOT.h"
#include <iostream>

#include "../../interface/GlobalVariables.h"
#include "../../interface/Samples/AllSamples.h"
#include "../../interface/Variables/Variable.h"
#include "../../interface/Objects/TdrStyle.h"

namespace std {

class multiFit{
public:
	multiFit();
	virtual ~multiFit();
	void allFits();
	void readHistos(AllSamples samples, Variable variable);
protected:
	TString objName;
	TString selection;
	TString folder;
	bool abseta, m3, mlb;
	void setSelection(TString sel_name);
	TLegend* legend(AllSamples samples);
	TText* doPrelim(double x_pos,double y_pos);
	TText* doChan(double x_pos,double y_pos);
	TH1D* readHistogram(Sample sample, Variable variable, bool btag);
	TH1D* readHistogram2(Sample sample, Variable variable, bool btag);
	THStack* buildStack(AllSamples samples, Variable variable);
	TH1D* allMChisto(AllSamples samples, Variable variable);
	TH1D* hashErrors(AllSamples samples, Variable variable);
	void addOverFlow(TH1D* overflow, Variable variable);
	void standardPlot(TH1D* data, THStack *hs, AllSamples samples, Variable variable);
	void ratioPlot(TH1D* data, THStack *hs, AllSamples samples, Variable variable);
	double readAndFit(AllSamples samples, Variable variable, TString syst_folder);
	double doFit(AllSamples samples, Variable variable, TString syst_folder);
	void drawTemplates(AllSamples samples, Variable variable, TString syst_folder);
	void normAndColor(TH1D* hist, Sample sample);
//	void fcn(int& npar, double* deriv, double& f, double par[], int flags);
};

class xSects{

	double AbsEta,M3,Mlb;
public : 
	xSects(double d, double e, double f){AbsEta=d; M3=e; Mlb=f;};
	double getAbsEta() {return AbsEta;};
	double getM3() {return M3;};
	double getMlb() {return Mlb;};
};

} /* namespace std */
#endif /* FIT_H_ */
