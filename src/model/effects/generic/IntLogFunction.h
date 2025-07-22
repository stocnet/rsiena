/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IntLogFunction.h
 *
 * Description: This file contains the definition of the
 * IntLogFunction class.
 *****************************************************************************/


#ifndef INTLOGFUNCTION_H_
#define INTLOGFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class LogTable;


class IntLogFunction: public AlterFunction
{
public:
	IntLogFunction(AlterFunction * pFunction);
	virtual ~IntLogFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter) const;

private:
	AlterFunction * lpFunction;
	LogTable * lpLogTable;
};

}

#endif /* INTLOGFUNCTION_H_ */
