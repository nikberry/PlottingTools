/*
 * GlobalVariables.h
 *
 *  Created on: Oct 24, 2013
 *      Author: philip
 */

#include "TString.h"

#ifndef GLOBALVARIABLES_H_
#define GLOBALVARIABLES_H_

namespace std {

struct Globals {
	static bool addRatioPlot;
	static bool addHashErrors;
	static bool addOverFlow;
	static TString lumi;
};

} /* namespace std */
#endif /* GLOBALVARIABLES_H_ */
