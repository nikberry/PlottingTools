/*
 * AllSamples.h
 *
 *  Created on: Oct 15, 2013
 *      Author: philip
 */

#ifndef ALLSAMPLES_H_
#define ALLSAMPLES_H_

#include "Sample.h"
#include "TObject.h"

namespace std {

class AllSamples {
public:
	AllSamples(TString systematic, TString eSystematic);
	virtual ~AllSamples();

	Sample* single_mu_data;
	Sample* ttbar;
	Sample* single_t;
	Sample* signal;
	Sample* vjets;
	Sample* qcd;
};

} /* namespace std */
#endif /* ALLSAMPLES_H_ */
