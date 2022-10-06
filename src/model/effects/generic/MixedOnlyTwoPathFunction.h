/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedOnlyTwoPath.h
 *
 * Description: This file contains the definition of the
 * ThreeCyclesFunction class.
 *****************************************************************************/

#ifndef MIXEDONLYTWOPATHFUNCTION_H_
#define MIXEDONLYTWOPATHFUNCTION_H_

#include "MixedNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;
// ----------------------------------------------------------------------------
// Section: MixedOnlyTwoPath class
// ----------------------------------------------------------------------------

/**
 * For sharedTo effect (see manual).
 */
class MixedOnlyTwoPath : public MixedNetworkAlterFunction, IntAlterFunction
{
public:
	MixedOnlyTwoPath(std::string dependentNetworkName,
			std::string explanatoryNetworkName, double parameter);

	virtual void initialize(const Data * pData,
			State * pState, int period, Cache * pCache);

	virtual double value(int alter);
	virtual int intValue(int alter);
};

}

#endif /* MIXEDONLYTWOPATHFUNCTION_H_ */
