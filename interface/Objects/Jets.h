/*
 * Jets.h
 *
 *  Created on: Oct 23, 2013
 *      Author: philip
 */

#include "../../interface/Objects/Object.h"

#ifndef JETS_H_
#define JETS_H_

namespace std {

class Jets: public Object {
public:
	Jets();
	virtual ~Jets();
	void allPlots(AllSamples samples);
};

} /* namespace std */
#endif /* JETS_H_ */
