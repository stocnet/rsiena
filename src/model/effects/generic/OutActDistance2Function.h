/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutActDistance2Function.h
 *
 * Description: This file contains the definition of the
 * OutActDistance2Function class.
 *****************************************************************************/

#ifndef OUTACTDISTANCE2FUNCTION_H_
#define OUTACTDISTANCE2FUNCTION_H_

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

class OutActDistance2Function: public MixedNetworkAlterFunction
{
public:
	OutActDistance2Function(std::string firstNetworkName,
			std::string secondNetworkName, double parameter, bool firstIn,
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
	std::string lvariableName {}; // local secondNetworkName
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
};

}

#endif /* OUTACTDISTANCE2FUNCTION_H_ */
