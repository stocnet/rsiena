/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ContinuousEffect.h
 *
 * Description: This file contains the definition of the
 * ContinuousEffect class.
 *****************************************************************************/

#ifndef CONTINUOUSEFFECT_H_
#define CONTINUOUSEFFECT_H_

#include "Effect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ContinuousLongitudinalData;


// ----------------------------------------------------------------------------
// Section: BehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects.
 */
class ContinuousEffect : public Effect
{
public:
	ContinuousEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual void preprocessEgo(int ego); // not used yet, only in internal functions

	/**
	 * Calculates how much this effect contributes to the change in the
	 * continuous behavior.
	 */
	virtual double calculateChangeContribution(int actor) = 0;

	inline double coefficient() const;
	void coefficient(double value);
	
	
	virtual double evaluationStatistic(double * currentValues); // not used yet
//	virtual double endowmentStatistic(const int * difference,
//		double * currentValues);
//	virtual double creationStatistic(int * difference,
//		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues); // not used yet
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues); // not used yet

protected:
	int n() const;
	double value(int actor) const;
	double centeredValue(int actor) const;
	bool missing(int observation, int actor) const;
	double range() const;
	double similarity(double a, double b) const;
	double similarityMean() const;

private:
	ContinuousLongitudinalData * lpContinuousData;
	const double * lvalues {};
	int lego {};
	
	// coefficient computed in the Bergstrom step by SdeSimulation, derived
	// from the parameter value
	double lcoefficient {};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the weight of this effect in the owner function, which is a weighted
 * sum of effects.
 */
double ContinuousEffect::coefficient() const
{
	return this->lcoefficient;
}

}

#endif /*CONTINUOUSEFFECT_H_*/
