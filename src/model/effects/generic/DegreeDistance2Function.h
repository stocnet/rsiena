/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DegreeDistance2Function.h
 *
 * Description: This file contains the definition of the
 * DegreeDistance2Function class.
 *****************************************************************************/

#ifndef DEGREEDISTANCE2FUNCTION_H_
#define DEGREEDISTANCE2FUNCTION_H_

#include "NetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

class DegreeDistance2Function: public NetworkAlterFunction
{
public:
	DegreeDistance2Function(std::string networkName,
			double parameter, bool firstIn,
			bool secondIn, bool average);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter) const;

private:
	bool lroot {}; // should the square root be taken?
	bool lfirstin {}; // first tie in- or outgoing?
	bool lsecondin {}; // second tie in- or outgoing?
	bool laverage {}; // should the average be used?
	double lavdegree {}; // average degree, all observations, secondNetworkName
	std::string lvariableName {}; // local networkName
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /* DEGREEDISTANCE2FUNCTION_H_ */
