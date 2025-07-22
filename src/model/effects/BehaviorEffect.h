/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * BehaviorEffect class.
 *****************************************************************************/

#ifndef BEHAVIOREFFECT_H_
#define BEHAVIOREFFECT_H_

#include "Effect.h"
#include <utility>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorLongitudinalData;


// ----------------------------------------------------------------------------
// Section: BehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects.
 */
class BehaviorEffect : public Effect
{
public:
	BehaviorEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData, State * pState,
			int period, Cache * pCache);
	virtual void initialize(const Data *pData, State *pState,
			State *pSimulatedState, int period, Cache *pCache);

	virtual void preprocessEgo(int ego);

	/**
	 * Calculates the change in the statistic corresponding to this effect if
	 * the given actor would change his behavior by the given amount.
	 */
	virtual double calculateChangeContribution(int actor,
			int difference) = 0;

	virtual double evaluationStatistic(double * currentValues);
	virtual std::pair<double,double *> evaluationStatistic(double * currentValues, bool needActorStatistics);
	virtual double endowmentStatistic(const int * difference,
			double * currentValues);
	virtual std::pair<double,double *> endowmentStatistic(const int * difference,
			double * currentValues, bool needActorStatistics);
	virtual double creationStatistic(int * difference,
			double * currentValues);
	virtual std::pair<double,double *> creationStatistic(int * difference,
			double * currentValues, bool needActorStatistics);
	virtual double egoStatistic(int ego, double * currentValues);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
			double * currentValues);

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
	virtual void initializeStatisticCalculation();
	virtual void cleanupStatisticCalculation();
	const int * initialValues();

private:
	BehaviorLongitudinalData * lpBehaviorData;
	const int * lvalues {};
	const int * linitialValues {};
	int lego {};
};

}

#endif /*BEHAVIOREFFECT_H_*/
