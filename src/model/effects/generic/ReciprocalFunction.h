/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocalFunction.h
 *
 * Description: This file contains the definition of the
 * ReciprocalFunction class.
 *****************************************************************************/

#ifndef RECIPROCALFUNCTION_H_
#define RECIPROCALFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class ReciprocalFunction: public AlterFunction
{
public:
	ReciprocalFunction(AlterFunction * pFunction);
	virtual ~ReciprocalFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter) const;

private:
	AlterFunction * lpFunction;
};

}

#endif /* RECIPROCALFUNCTION_H_ */
