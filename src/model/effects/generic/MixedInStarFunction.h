/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedInStarFunction.h
 *
 * Description: This file contains the definition of the
 * MixedInStarFunction class.
 *****************************************************************************/

#ifndef MIXEDINSTARFUNCTION_H_
#define MIXEDINSTARFUNCTION_H_

#include "MixedNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

class MixedConfigurationTable;


/**
 * Defines a function that returns the number of two-paths between
 * the ego and alters in a network of the given name.
 */
class MixedInStarFunction:
public MixedNetworkAlterFunction, IntAlterFunction
{
public:
	MixedInStarFunction(std::string firstNetworkName, std::string secondNetworkName);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter);
	virtual int intValue(int alter);

private:
	MixedConfigurationTable * lpTable;
};

}

#endif /* MIXEDINSTARFUNCTION_H_ */
