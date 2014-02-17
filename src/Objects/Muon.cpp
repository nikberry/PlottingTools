/*
 * Muon.cpp
 *
 *  Created on: Dec 24, 2013
 *      Author: philip
 */

#include "../../interface/Objects/Muon.h"
#include "TCanvas.h"

namespace std {

Muon::Muon() {
	objName = "Muon";

}

Muon::~Muon() {

}

void Muon::allPlots(AllSamples samples){

	Variable pt("muon_pT", "muon p_{T}", 0, 300, 8); //2 for all default
	savePlot(samples, pt);

	Variable absEta("muon_AbsEta", "muon |#eta|", 0, 2.6, 8);
	savePlot(samples, absEta);

	Variable eta("muon_eta", "muon #eta", -2.6, 2.6, 10);
	savePlot(samples, eta);

	Variable phi("muon_phi", "muon #phi", -3.4, 3.4, 8);
	savePlot(samples, phi);

}

} /* namespace std */
