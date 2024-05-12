/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreeActivityEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreeActivityEffect class.
 *****************************************************************************/

#ifndef OUTDEGREEACTIVITYEFFECT_H_
#define OUTDEGREEACTIVITYEFFECT_H_

#include <string>
#include "NetworkEffect.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * This class defines the outdegree activity effect defined by
 * s_i(x)= x_{i+}^2. The corresponding statistic is
 * the sum of squared outdegrees over all actors.
 */
class OutdegreeActivityEffect : public NetworkEffect
{
friend class BothDegreesEffect;

public:
	OutdegreeActivityEffect(const EffectInfo * pEffectInfo, bool centered);

	virtual void initialize(const Data * pData, State * pState,	int period,
			Cache * pCache);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	virtual double endowmentStatistic(Network * pLostTieNetwork);

private:
	bool lcentered {};
	double lcentering {};
	string lvariableName {};

};

}

#endif /*OUTDEGREEACTIVITYEFFECT_H_*/
