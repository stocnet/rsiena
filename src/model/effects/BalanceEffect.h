/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BalanceEffect.h
 *
 * Description:  This file contains the declaration of the class
 * BalanceEffect.
 *****************************************************************************/

#ifndef BALANCEEFFECT_H_
#define BALANCEEFFECT_H_

#include "NetworkEffect.h"
#include "network/IncidentTieIterator.h"

namespace siena
{

/**
 * This class defines the balance effect.
 */
class BalanceEffect : public NetworkEffect
{
public:
	BalanceEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double calculateContribution(int alter) const;

protected:
	virtual void initializeStatisticCalculation();
	virtual void cleanupStatisticCalculation();
	virtual double tieStatistic(int alter);

private:
	void markInvalidActors(IncidentTieIterator iter,
		int & validActorCount);

	// The centering constant b_0.
	double lbalanceMean {};

	// An indicator array for invalid actors (used in statistic calculations)
	// Invariants:
	// A: lflag[i] <= lround for all actors
	// B: lflag[i] == lround for invalid actors

	int * lflag {};
	int lround {};
};

}

#endif /*BALANCEEFFECT_H_*/
