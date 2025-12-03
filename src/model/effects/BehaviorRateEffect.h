/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorRateEffect.h
 *
 * Description: This file contains the definition of the
 * BehaviorRateEffect class.
 *****************************************************************************/

#ifndef BEHAVIORERATEFFECT_H_
#define BEHAVIORERATEFFECT_H_

#include "Effect.h"
#include <utility>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorLongitudinalData;


// ----------------------------------------------------------------------------
// Section: BehaviorRateEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects.
 */
class BehaviorRateEffect : public Effect
{
public:
	BehaviorRateEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData, State * pState,
			int period, Cache * pCache);

	virtual double calculateContribution(int i) const;


protected:
	int n() const;
	int value(int actor) const;
	int initialValue(int actor) const;
	double centeredValue(int actor) const;
	double overallCenterMean() const;
	bool missing(int observation, int actor) const;
	double range() const;
	double similarity(double a, double b) const;
	double similarityMean() const;
    double variance() const;
	const int * initialValues();

private:
	BehaviorLongitudinalData * lpBehaviorData;
	const int * lvalues {};
	const int * linitialValues {};
};

}

#endif /*BEHAVIORRATEFFECT_H_*/
