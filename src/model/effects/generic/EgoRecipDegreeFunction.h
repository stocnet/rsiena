/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoRecipDegreeFunction.h
 *
 * Description: This file contains the definition of the
 * EgoRecipDegreeFunction class.
 *****************************************************************************/

#ifndef EGORECIPDEGREEFUNCTION_H_
#define EGORECIPDEGREEFUNCTION_H_

#include "OneModeNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the reciprocal degree of the ego regardless
 * of the alter.
 */
class EgoRecipDegreeFunction: public OneModeNetworkAlterFunction, IntAlterFunction
{
public:
	EgoRecipDegreeFunction(std::string networkName); 
	
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter) const;
	virtual int intValue(int alter);
};

}

#endif /* EGORECIPDEGREEFUNCTION_H_ */
