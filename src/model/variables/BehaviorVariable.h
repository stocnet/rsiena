/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorVariable.h
 *
 * Description: This file contains the definition of the
 * BehaviorVariable class.
 *****************************************************************************/

#ifndef BEHAVIORVARIABLE_H_
#define BEHAVIORVARIABLE_H_

#include "DependentVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorLongitudinalData;
class SimulationActorSet;


// ----------------------------------------------------------------------------
// Section: BehaviorVariable class
// ----------------------------------------------------------------------------

/**
 * This class represents the state of a behavioral dependent variable.
 * @see DependentVariable
 */
class BehaviorVariable : public DependentVariable
{
public:
	BehaviorVariable(BehaviorLongitudinalData * pData,
		EpochSimulation * pSimulation);
	virtual ~BehaviorVariable();

	virtual int m() const;
	virtual LongitudinalData * pData() const;
	virtual bool behaviorVariable() const;
	virtual void initialize(int period);
	virtual void setLeaverBack(const SimulationActorSet * pActorSet,
		int actor);

	void behaviorModelType(int type);
	virtual BehaviorModelType behaviorModelType() const;
	virtual bool behaviorModelTypeABSORB() const;

	virtual void makeChange(int actor);

	bool missingStartValue(int actor) const;
	int value(int actor) const;
	void value(int actor, int newValue);
	double centeredValue(int actor) const;
	bool structural(int actor) const;
	double similarity(int i, int j) const;
	const int * values() const;
	int range() const;
	double similarityMean() const;

	virtual double probability(MiniStep * pMiniStep);
	virtual bool validMiniStep(const MiniStep * pMiniStep,
		bool checkUpOnlyDownOnlyConditions = true) const;
	virtual MiniStep * randomMiniStep(int ego);
	virtual bool missing(const MiniStep * pMiniStep) const;
	virtual bool structural(const MiniStep * pMiniStep) const;

private:
	void preprocessEgo();
	void preprocessEffects(const Function * pFunction);
	double totalEvaluationContribution(int actor,
		int difference) const;
	double totalEndowmentContribution(int actor,
		int difference) const;
	double totalCreationContribution(int actor,
		int difference) const;
	void accumulateScores(int difference, bool UpPossible,
		bool downPossible) const;
	void calculateProbabilities(int actor);
	void accumulateDerivatives() const;

	// The observed data for this behavioral variable
	BehaviorLongitudinalData * lpData;

	// The current value of the variable per each actor
	int * lvalues {};

	// A two-dimensional array of change contributions to effects, where
	// rows correspond to differences and columns correspond to effects in the
	// evaluation function.
	double ** levaluationEffectContribution {};

	// A two-dimensional array of change contributions to effects, where
	// rows correspond to differences and columns correspond to effects in the
	// endowment function.
	double ** lendowmentEffectContribution {};

	// A two-dimensional array of change contributions to effects, where
	// rows correspond to differences and columns correspond to effects in the
	// tie creation function.
	double ** lcreationEffectContribution {};

	// Selection probability per each difference:
	// lprobabilities[0] - probability for a downward change
	// lprobabilities[1] - probability of no change
	// lprobabilities[2] - probability of an upward change
	double * lprobabilities {};

	// Indicates if upward change is possible in the current situation
	bool lupPossible {};

	// Indicates if downward change is possible in the current situation
	bool ldownPossible {};

	// the actor under consideration
	int lego {};

	// the model type
	BehaviorModelType lbehaviorModelType {};
};

}

#endif /*BEHAVIORVARIABLE_H_*/
