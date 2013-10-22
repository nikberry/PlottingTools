/*
 * Variable.h
 *
 *  Created on: Oct 17, 2013
 *      Author: philip
 */

#include "TString.h"

#ifndef VARIABLE_H_
#define VARIABLE_H_

namespace std {

class Variable {
public:
	Variable();
	Variable(TString name_temp, TString xTitle_temp, double minX, double maxX, int rebinFact_temp);
	virtual ~Variable();
	TString name;
	TString xTitle;
	double minX;
	double maxX;
	int rebinFact;
};

} /* namespace std */
#endif /* VARIABLE_H_ */
