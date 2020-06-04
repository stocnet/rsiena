/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DifferenceFunction.h
 *
 * Description: This file contains the definition of the
 * DifferenceFunction class.
 *****************************************************************************/

#ifndef DIFFERENCEFUNCTION_H_
#define DIFFERENCEFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class DifferenceFunction: public AlterFunction
{
public:
	DifferenceFunction(AlterFunction * pFirstFunction,
		AlterFunction * pSecondFunction);
	virtual ~DifferenceFunction();

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

#endif /* DIFFERENCEFUNCTION_H_ */
