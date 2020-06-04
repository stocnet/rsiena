/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoFunction.h
 *
 * Description: This file contains the definition of the
 * EgoFunction class.
 *****************************************************************************/

#ifndef EGOFUNCTION_H_
#define EGOFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class EgoFunction: public AlterFunction
{
public:
	EgoFunction(AlterFunction * pFirstFunction);
	virtual ~EgoFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter);

private:
	AlterFunction * lpFirstFunction;
};

}

#endif /* EGOFUNCTION_H_ */
