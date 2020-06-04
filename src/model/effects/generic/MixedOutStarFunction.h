/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedOutStarFunction.h
 *
 * Description: This file contains the definition of the
 * MixedOutStarFunction class.
 *****************************************************************************/

#ifndef MIXEDOUTSTARFUNCTION_H_
#define MIXEDOUTSTARFUNCTION_H_

#include "MixedNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

class MixedConfigurationTable;


/**
 * Defines a function that returns the number of two-paths between
 * the ego and alters in a network of the given name.
 */
class MixedOutStarFunction:
public MixedNetworkAlterFunction, IntAlterFunction
{
public:
	MixedOutStarFunction(std::string firstNetworkName,
			std::string secondNetworkName);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter);
	virtual int intValue(int alter);

private:
	MixedConfigurationTable * lpTable;
};

}

#endif /* MIXEDOUTSTARFUNCTION_H_ */
