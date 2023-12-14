/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DenseTriadsSimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeEffect class.
 *****************************************************************************/

#ifndef DENSETRIADSSIMILARITYEFFECT_H_
#define DENSETRIADSSIMILARITYEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Similarity in dense triads effect defined as
 * s_i(x) = sum_{j,h} (sim(v_i, v_j) + sim(v_i, v_h)) x
 *   I{x_{ij} + x_{ji} + x_{ih} + x_{hi} + x_{jh} + x_{hj} >= c},
 * where sim(v_i, v_j) are the centered similarity scores.
 */
class DenseTriadsSimilarityEffect : public NetworkDependentBehaviorEffect
{
public:
	DenseTriadsSimilarityEffect(const EffectInfo * pEffectInfo);
	virtual ~DenseTriadsSimilarityEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);

private:
	void updateMarks(int i);

	int ldensity {};

	// A helper array of marks
	int * lmark {};

	// Given an ego i
	// mark[h] = baseMark + 2 if there are mutual ties between i and h,
	// mark[h] = baseMark + 1 if only one of the mutual ties is present,
	// mark[h] <= baseMark otherwise.

	int lbaseMark {};
};

}

#endif /*DENSETRIADSSIMILARITYEFFECT_H_*/

