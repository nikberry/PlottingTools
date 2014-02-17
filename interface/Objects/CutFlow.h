/*
 * CutFlow.h
 *
 *  Created on: Oct 26, 2013
 *      Author: philip
 */

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"

#include <iostream>

#include "../../interface/GlobalVariables.h"
#include "../../interface/Samples/AllSamples.h"
#include "../../interface/Variables/Variable.h"
#include "../../interface/Objects/TdrStyle.h"

#ifndef CUTFLOW_H_
#define CUTFLOW_H_

namespace std {

class CutFlow {
public:
	CutFlow();
	virtual ~CutFlow();
	void allPlots(AllSamples samples);
private:
	TString objName;
	TString selection;
	void setSelection(TString sel_name);
	TLegend* legend(AllSamples samples);
	TText* doPrelim(double x_pos,double y_pos);
	TText* doChan(double x_pos,double y_pos);
	THStack* buildStack(AllSamples samples, Variable variable);
	TH1D* allMChisto(AllSamples samples, Variable variable);
	TH1D* hashErrors(AllSamples samples, Variable variable);
	TH1D* readCutFlowHistogram(Sample sample, Variable variable);
	void readCutFlowHistos(AllSamples samples, Variable variable);
	void saveCutFlowPlot(AllSamples samples, Variable variable);
	void setBinLabels(THStack* hs, TH1D* data);
	void cutFlowEffPlots(TH1D* data, AllSamples samples, Variable variable);
	void standardCutFlowPlot(TH1D* data, THStack *hs, AllSamples samples, Variable variable);
	void ratioCutFlowPlot(TH1D* data, THStack *hs, AllSamples samples, Variable variable);
	void writeTable(AllSamples samples, Variable variable);
};

} /* namespace std */
#endif /* CUTFLOW_H_ */
