/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SumFunction.h
 *
 * Description: This file contains the definition of the
 * SumFunction class.
 *****************************************************************************/

#ifndef SUMFUNCTION_H_
#define SUMFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class SumFunction: public AlterFunction
{
public:
	SumFunction(AlterFunction * pFirstFunction,
		AlterFunction * pSecondFunction);
	virtual ~SumFunction();

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

#endif /* SUMFUNCTION_H_ */
