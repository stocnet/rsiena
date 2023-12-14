/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DenseTriadsBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeEffect class.
 *****************************************************************************/

#ifndef DENSETRIADSBEHAVIOREFFECT_H_
#define DENSETRIADSBEHAVIOREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Dense triads behavior effect (see manual).
 */
class DenseTriadsBehaviorEffect : public NetworkDependentBehaviorEffect
{
public:
	DenseTriadsBehaviorEffect(const EffectInfo * pEffectInfo);
	virtual ~DenseTriadsBehaviorEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	int denseTriadCount(int i);

	int ldensity {};

	// A helper array of marks for counting dense triads
	int * lmark {};

	// Given an ego i
	// mark[h] = baseMark + 2 if there are mutual ties between i and h,
	// mark[h] = baseMark + 1 if only one of the mutual ties is present,
	// mark[h] <= baseMark otherwise.

	int lbaseMark {};
};

}

#endif /*DENSETRIADSBEHAVIOREFFECT_H_*/

