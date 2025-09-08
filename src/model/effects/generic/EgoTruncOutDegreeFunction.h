/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File:EgoTruncOutDegreeFunction.h
 *
 * Description: This file contains the definition of the
 *EgoTruncOutDegreeFunction class.
 *****************************************************************************/

#ifndef EGOTRUNCOUTDEGREEFUNCTION_H_
#define EGOTRUNCOUTDEGREEFUNCTION_H_


#include <string>
#include "NetworkAlterFunction.h" 


namespace siena
{
// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: EgoTruncOutDegreeFunction class
// ----------------------------------------------------------------------------


/**
 * Defines a function that returns the out-degree of the ego regardless
 * of the alter.
 */
class EgoTruncOutDegreeFunction: public NetworkAlterFunction
{
public:
	EgoTruncOutDegreeFunction(std::string networkName, bool root, 
								bool threshold, int p);
	virtual double value(int alter) const;
	
private:
	// Indicates if the square root of outdegrees must be used
	bool lroot {};
	bool lthreshold {};

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
	double lp {};
	int lintp {};
	std::string lvariableName {};
};

}

#endif /*EGOTRUNCOUTDEGREEFUNCTION_H_ */
