/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EqualsZeroFunction.h
 *
 * Description: This file contains the definition of the
 * EqualsZeroFunction class.
 * It is currently not used, and therefore not included in sources.list,
 *****************************************************************************/


#ifndef EQUALSZEROFUNCTION_H_
#define EQUALSZEROFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{


class EqualsZeroFunction: public AlterFunction
{
public:
	EqualsZeroFunction(AlterFunction * pFunction, bool minus);
	virtual ~EqualsZeroFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter) const;

private:
	AlterFunction * lpFunction;
	bool lminus {};
};

}

#endif /* EQUALSZEROFUNCTION_H_ */
