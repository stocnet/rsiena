/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedTwoPathFunction.h
 *
 * Description: This file contains the definition of the
 * MixedTwoPathFunction class.
 *****************************************************************************/

#ifndef MIXEDTWOSTEPFUNCTION_H_
#define MIXEDTWOSTEPFUNCTION_H_

#include "MixedNetworkAlterFunction.h"
#include "IntAlterFunction.h"
#include "network/NetworkUtils.h"

namespace siena
{

class MixedConfigurationTable;


/**
 * Defines a function that returns the number of two-paths between
 * the ego and alters in a network of the given name.
 */
class MixedTwoStepFunction:
public MixedNetworkAlterFunction, IntAlterFunction
{
public:
	MixedTwoStepFunction(string firstNetworkName, string secondNetworkName,
						Direction firstDirection, Direction secondDirection);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter);
	virtual int intValue(int alter);

private:
	MixedConfigurationTable * lpTable;

	Direction ldirection1;
	Direction ldirection2;
};

}

#endif /* MIXEDTWOSTEPFUNCTION_H_ */
