/*
 * Object.h
 *
 *  Created on: Oct 18, 2013
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

#include "../../interface/Samples/AllSamples.h"
#include "../../interface/Variables/Variable.h"
#include "../../interface/Objects/TdrStyle.h"

#ifndef OBJECT_H_
#define OBJECT_H_

namespace std {

class Object {
public:
	Object();
	virtual ~Object();
	void allPlots(AllSamples samples);
protected:
	TLegend* legend(AllSamples samples);
	TText* doPrelim(double x_pos,double y_pos);
	TText* doChan(double x_pos,double y_pos);
	TH1D* readGe2bHistogram(Sample sample, Variable variable);
	THStack* buildStack(AllSamples samples, Variable variable);
	void savePlot(AllSamples samples, Variable variable);
};

} /* namespace std */
#endif /* OBJECT_H_ */
