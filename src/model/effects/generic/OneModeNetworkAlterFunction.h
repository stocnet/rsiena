/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OneModeNetworkAlterFunction.h
 *
 * Description: This file contains the definition of the
 * OneModeNetworkAlterFunction class.
 *****************************************************************************/

#ifndef ONEMODENETWORKALTERFUNCTION_H_
#define ONEMODENETWORKALTERFUNCTION_H_

#include "NetworkAlterFunction.h"

namespace siena
{

/**
 * An alter function that depends on a one-mode network.
 */
class OneModeNetworkAlterFunction: public NetworkAlterFunction
{
public:
	OneModeNetworkAlterFunction(std::string networkName);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
};

}

#endif /* ONEMODENETWORKALTERFUNCTION_H_ */
