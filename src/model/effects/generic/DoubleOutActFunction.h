/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleOutActFunction.h
 *
 * Description: This file contains the definition of the
 * ThreeCyclesFunction class.
 *****************************************************************************/

#ifndef DOUBLEOUTACTFUNCTION_H_
#define DOUBLEOUTACTFUNCTION_H_

#include "MixedNetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;
// ----------------------------------------------------------------------------
// Section: DoubleOutActFunction class
// ----------------------------------------------------------------------------

/**
 * For sharedTo effect (see manual).
 */
class DoubleOutActFunction : public MixedNetworkAlterFunction
{
public:
	DoubleOutActFunction(std::string dependentNetworkName,
			std::string explanatoryNetworkName, double parameter, bool change);

	virtual void initialize(const Data * pData,
			State * pState, int period, Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lroot {}; // should the square root be taken?
	bool lchange {}; // should the change statistic be calculated?
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /*DOUBLEOUTACTFUNCTION_H_*/
