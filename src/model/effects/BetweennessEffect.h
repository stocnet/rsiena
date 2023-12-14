/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BetweennessEffect.h
 *
 * Description: This file contains the declaration of the class
 * BetweennessEffect.
 *****************************************************************************/

#ifndef BETWEENNESSEFFECT_H_
#define BETWEENNESSEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the betweenness effect.
 */
class BetweennessEffect : public NetworkEffect
{
public:
	BetweennessEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual void initializeStatisticCalculation();
	virtual void onNextEgo(int ego);
	virtual double tieStatistic(int alter);
	virtual void cleanupStatisticCalculation();

private:
	// Helper array of marks for statistic calculations

	int * lmark {};
	int lcurrentMark {};
	int lbaseMark {};
};

}

#endif /*BETWEENNESSEFFECT_H_*/
