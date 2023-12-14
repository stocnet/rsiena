/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AntiIsolateEffect.h
 *
 * Description: This file contains the definition of the
 * AntiIsolateEffect class.
 *****************************************************************************/

#ifndef AntiIsolateEffect_H_
#define AntiIsolateEffect_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the anti in-isolates effect defined by
 * the change statistic -1 if x_{ij} = x_{+j} = 1
 * and +1 if x_{ij} = x_{+j} = 0.
 * This is not an evaluation/endowment/creation type of effect,
 * but a gratification effect.
 * The target statistic is the number of actors with indegree at least 1.
 */
class AntiIsolateEffect : public NetworkEffect
{
public:
	AntiIsolateEffect(const EffectInfo * pEffectInfo,
		bool outAlso,
		int minDegree);
	virtual double calculateContribution(int alter) const;
// The minDegree could also be dealt with by an effect parameter.
// But I (t.s.) wish to include in any case the minDegree 1 and 2 and 3 cases
// as separate effects, and using it for minDegree 4 and more seems overdone.
// Although, at least 3 is a group, so 3 might also be relevant.

protected:
	virtual double egoStatistic(int ego, const Network * pNetwork);

private:
	// Indicates the minimum degree threshold
	int lminDegree {};
	// Indicates that also out-tie-isolation of alter is considered
	bool loutAlso {};
};

}

#endif /*AntiIsolateEffect_H_*/
