/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IntSqrtFunction.h
 *
 * Description: This file contains the definition of the
 * IntSqrtFunction class.
 *****************************************************************************/


#ifndef INTSQRTFUNCTION_H_
#define INTSQRTFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class SqrtTable;


class IntSqrtFunction: public AlterFunction
{
public:
	IntSqrtFunction(AlterFunction * pFunction);
	virtual ~IntSqrtFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter);

private:
	AlterFunction * lpFunction;
	SqrtTable * lpSqrtTable;
};

}

#endif /* INTSQRTFUNCTION_H_ */
