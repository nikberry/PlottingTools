/*
 * MET.cpp
 *
 *  Created on: Oct 11, 2013
 *      Author: philip
 */

#include "../../interface/Objects/MET.h"
#include "TCanvas.h"

namespace std {

MET::MET() {
	objName = "MET";
}

MET::~MET() {

}

void MET::allPlots(AllSamples samples){

	setMetType("patType1CorrectedPFMet");

//  Variable::Variable(TString name_temp, TString xTitle_temp, double minX_temp, double maxX_temp, int rebinFact_temp)
	Variable met;
	savePlot(samples, met);

	Variable met_sig("METsignificance", "E_{T}^{miss} significance", 0, 150, 2);
	savePlot(samples, met_sig);

	Variable st("ST", "ST [GeV]", 0, 1500, 2);
	savePlot(samples, st);

	Variable wpt("WPT", "p_{T}(W) [GeV]", 0, 300, 5);
	savePlot(samples, wpt);

	Variable mt("MT", "m_{T}(W) [GeV]", 0, 300, 5);
	savePlot(samples, mt);

	Variable met_phi("MET_phi", "E_{T}^{miss} #phi", -3.4, 3.4, 1);
	savePlot(samples, met_phi);

}

void MET::setMetType(TString metType){
	objName += "/"+metType;
}

} /* namespace std */


