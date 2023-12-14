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

#ifndef CONTINUOUSVARIABLE_H_
#define CONTINUOUSVARIABLE_H_

#include <vector>
#include "model/Function.h"
#include "utils/NamedObject.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ContinuousLongitudinalData;
class EpochSimulation;
class SimulationActorSet;


// ----------------------------------------------------------------------------
// Section: BehaviorVariable class
// ----------------------------------------------------------------------------

/**
 * This class is a DUMMY state of a continuous behavioral dependent variable.
 * @see DependentVariable
 */
class ContinuousVariable : public NamedObject
{
public:
	ContinuousVariable(ContinuousLongitudinalData * pData,
		EpochSimulation * pSimulation);
	virtual ~ContinuousVariable();
	
	// Initialization
	void initialize(int period);
	void initializeFunction() const;
		
	// Accessors
	int n() const;
	int id() const;	// not yet used
	double simulatedDistance() const; // not yet used
	double value(int actor) const; // not yet used
	
	inline const double * values() const;
	inline const Function * pFunction() const;

	// Computation
	void calculateEffectContribution();
	double totalFunctionContribution(int actor) const;
	void accumulateScores(const vector<double> &actorMeans, 
	                      const vector<double> &actorErrors,
						  double dt);
	
	// Setters
	void value(int actor, double newValue);

protected:
	inline EpochSimulation * pSimulation() const;
	void simulatedDistance(double distance); // not yet used
	
private:
	// A simulation of the actor-based model, which owns this variable
	EpochSimulation * lpSimulation;

	// The underlying set of actors
	const SimulationActorSet * lpActorSet;

	// The current period (in [0, observations - 2])
	int lperiod {};

	// The observed data for this continuous behavioral variable
	ContinuousLongitudinalData * lpData;

	// The current value of the variable per each actor
	double * lvalues {};

	// The distance of this variable from the observed data at the beginning
	// of the current period
	double lsimulatedDistance {};

	// A storage box for all the effect in the model for this variable
	Function * lpFunction;
	
	// The basic scale parameter for the current period
	double lbasicScale {};

	// The score for the basic scale parameter for this variable for this period
	double lbasicScaleScore {};

	// The derivative for the basic rate parameter for this variable for
	// this period
	double lbasicScaleDerivative {};
	
	// A two-dimensional array of contributions of effects, where rows correspond 
	// to actors and columns to effects in the SDE (needed for determining the
	// score function)
	double ** leffectContribution {};

};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the actor-based model owning this variable.
 */
EpochSimulation * ContinuousVariable::pSimulation() const
{
	return this->lpSimulation;
}

/**
 * Returns the array of current values for this variable.
 */
const double * ContinuousVariable::values() const
{
	return this->lvalues;
}

/**
 * Returns the storage box with all effects for this variable.
 */
const Function * ContinuousVariable::pFunction() const
{
	return this->lpFunction;
}

}

#endif /*CONTINUOUSVARIABLE_H_*/
