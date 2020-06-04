/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConditionalFunction.h
 *
 * Description: This file contains the definition of the
 * ConditionalFunction class.
 *****************************************************************************/

#ifndef CONDITIONALFUNCTION_H_
#define CONDITIONALFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class AlterPredicate;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines an alter function depending on a condition. More specifically,
 * this function owns two sub-functions and a predicate. If the predicate
 * is true then the calculation of the function value is delegated to the
 * if-function. Otherwise, the calculation is delegated to the else-function.
 */
class ConditionalFunction: public AlterFunction
{
public:
	ConditionalFunction(AlterPredicate * pPredicate,
		AlterFunction * pIfFunction,
		AlterFunction * pElseFunction);
	virtual ~ConditionalFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter);

private:
	AlterPredicate * lpPredicate;
	AlterFunction * lpIfFunction;
	AlterFunction * lpElseFunction;
};

}

#endif /* CONDITIONALFUNCTION_H_ */
