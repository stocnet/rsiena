/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BothDegreesEffect.h
 *
 * Description: This file contains the definition of the
 * BothDegreesEffect class.
 *****************************************************************************/

#ifndef BOTHDEGREESEFFECT_H_
#define BOTHDEGREESEFFECT_H_

#include <string>
#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * This class defines the both degrees effect defined as
 * the sum of the indegree popularity and outdegree activity effects.
 */
class BothDegreesEffect : public NetworkEffect
{
public:
	BothDegreesEffect(const EffectInfo * pEffectInfo, bool centered);

	virtual void initialize(const Data * pData, State * pState,	int period,
			Cache * pCache);
	virtual double calculateContribution(int alter) const;
	virtual double endowmentStatistic(Network * pLostTieNetwork);


protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the square root of degrees must be used
	bool lroot {};
	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
	bool lcentered {};
	double lcentering {};
	std::string lvariableName {};
};

}

#endif /*BOTHDEGREESEFFECT_H_*/
