/*
 * Variable.cpp
 *
 *  Created on: Oct 17, 2013
 *      Author: philip
 */

#include "../../interface/Variables/Variable.h"

namespace std {

Variable::Variable() {
	name = "MET";
	xTitle = "E_{T}^{miss}";
	minX = 0.;
	maxX = 300;
	rebinFact = 2;
}

Variable::Variable(TString name_temp, TString xTitle_temp, double minX_temp, double maxX_temp, int rebinFact_temp) {
	name = name_temp;
	xTitle = xTitle_temp;
	minX = minX_temp;
	maxX = maxX_temp;
	rebinFact = rebinFact_temp;
}

Variable::~Variable() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
