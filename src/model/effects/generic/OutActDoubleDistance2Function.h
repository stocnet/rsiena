/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutActDistance2Function.h
 *
 * Description: This file contains the definition of the
 * OutActDoubleDistance2Function class.
 *****************************************************************************/

#ifndef OUTACTDOUBLEDISTANCE2FUNCTION_H_
#define OUTACTDOUBLEDISTANCE2FUNCTION_H_

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

class OutActDoubleDistance2Function: public MixedNetworkAlterFunction
{
public:
	OutActDoubleDistance2Function(std::string firstNetworkName,
			std::string secondNetworkName, double parameter,
			bool secondIn, bool average);
	virtual ~OutActDoubleDistance2Function();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
	double increase(int h) const;
	virtual void preprocessEgo(int ego);

	virtual double value(int alter) const;

private:
	bool lroot {}; // should the square root be taken?
	bool lsecondin {}; //  in- or outdegrees?
	bool laverage {}; // should the average be used?
	double lavdegree {}; // average degree, all observations, secondNetworkName
	std::string lvariableName {}; // local secondNetworkName
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
	// for use in preprocessing:
	int * ltimesFound {};
	double loutInDist2Degree {};
	int ln {};
};

}

#endif /* OUTACTDOUBLEDISTANCE2FUNCTION_H_ */
