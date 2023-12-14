/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageDegreeEffect.h
 *
 * Description: This file contains the definition of the
 * AverageDegreeEffect class.
 *****************************************************************************/

#ifndef AVERAGEDEGREEEFFECT_H_
#define AVERAGEDEGREEEFFECT_H_

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
 * This class defines the average degree effect
   defined by the actor's outdegree multiplied by the average degree.
 */
class AverageDegreeEffect : public NetworkEffect
{
public:
	AverageDegreeEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData, State * pState,	int period,
			Cache * pCache);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego, const Network * pNetwork);
private:
	double lcentering {};
};

}

#endif /*AVERAGEDEGREEEFFECT_H_*/
