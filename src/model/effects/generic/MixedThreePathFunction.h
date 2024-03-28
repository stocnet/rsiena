/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedThreePathFunction.h
 *
 * Description: This file contains the definition of the
 * MixedThreePathFunction class.
 *****************************************************************************/

#ifndef MIXEDTHREEPATHFUNCTION_H_
#define MIXEDTHREEPATHFUNCTION_H_

#include "MixedNetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

class MixedThreePathFunction: public MixedNetworkAlterFunction
{
public:
	MixedThreePathFunction(std::string firstNetworkName,
			std::string secondNetworkName, double parameter,
			bool firstIn, bool secondIn, bool average);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lroot {}; // should the square root be taken?
	bool lcenter {}; // should there be centering?
	bool lfirstin {}; // first tie in- or outgoing?
	bool lsecondin {}; // second tie in- or outgoing?
	bool laverage {}; // should the average be taken?
	double lavdegree {}; // if centering: average degree, all observations, secondNetworkName
	std::string lvariableName {}; // local secondNetworkName
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /* MIXEDTHREEPATHFUNCTION_H_ */
