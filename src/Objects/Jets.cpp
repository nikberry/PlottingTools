/*
 * Jets.cpp
 *
 *  Created on: Oct 23, 2013
 *      Author: philip
 */

#include "../../interface/Objects/Jets.h"

namespace std {

Jets::Jets() {
	// TODO Auto-generated constructor stub
	objName = "Jets";
}

Jets::~Jets() {
	// TODO Auto-generated destructor stub
}

void Jets::allPlots(AllSamples samples){

//  Variable::Variable(TString name_temp, TString xTitle_temp, double minX_temp, double maxX_temp, int rebinFact_temp)
	Variable pt_all("all_jet_pT", "all jet p_{T}", 0, 300, 5);
	savePlot(samples, pt_all);

	Variable phi_all("all_jet_phi", "all jet #phi", -3.4, 3.4, 5);
	savePlot(samples, phi_all);

	Variable eta_all("all_jet_eta", "all jet #eta", -2.6, 2.6, 5);
	savePlot(samples, eta_all);
}

} /* namespace std */
