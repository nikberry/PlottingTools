/*
 * Muon.h
 *
 *  Created on: Dec 24, 2013
 *      Author: philip
 */

#include "../../interface/Objects/Object.h"

#ifndef MUON_H_
#define MUON_H_

namespace std {

class Muon: public Object {
public:
	Muon();
	virtual ~Muon();
	void allPlots(AllSamples samples);
};

} /* namespace std */
#endif /* MUON_H_ */
