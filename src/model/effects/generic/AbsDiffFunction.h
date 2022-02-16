/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DifferenceFunction.h
 *
 * Description: This file contains the definition of the
 * AbsDiffFunction class.
 *****************************************************************************/

#ifndef ABSDIFFFUNCTION_H_
#define ABSDIFFFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class AbsDiffFunction: public AlterFunction
{
public:
	AbsDiffFunction(AlterFunction * pFirstFunction,
		AlterFunction * pSecondFunction);
	virtual ~AbsDiffFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter);

private:
	AlterFunction * lpFirstFunction;
	AlterFunction * lpSecondFunction;
};

}

#endif /* ABSDIFFFUNCTION_H_ */
