/*
 * MET.h
 *
 *  Created on: Oct 11, 2013
 *      Author: philip
 */

#include "../../interface/Objects/Object.h"

#ifndef MET_H_
#define MET_H_

namespace std {

class MET: public Object {
public:
	MET();
	virtual ~MET();
	void allPlots(AllSamples samples);
private:
	void setMetType(TString metType);
};

} /* namespace std */
#endif /* MET_H_ */
