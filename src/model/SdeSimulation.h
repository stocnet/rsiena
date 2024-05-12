/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SdeSimulation.h
 *
 * Description: This file contains the definition of the
 * SdeSimulation class.
 *****************************************************************************/

#ifndef SDESIMULATION_H_
#define SDESIMULATION_H_

#include <vector>

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

//class ContinuousLongitudinalData;
class EpochSimulation;
//class MatrixXd;


// ----------------------------------------------------------------------------
// Section: SdeSimulation class
// ----------------------------------------------------------------------------

/**
 * This class is a DUMMY SdeSimulation
 * @see SdeSimulation
 */
class SdeSimulation
{
public:
	SdeSimulation(EpochSimulation * pSimulation);
	virtual ~SdeSimulation();

	inline double feedbackParameter() const;
	inline double feedbackCoefficient() const;
	inline double wienerParameter() const;
	inline double wienerCoefficient() const; 
	
	void initialize(int period);
	void setBergstromCoefficients(double dt);
	double randomComponent() const;
	
	double basicScaleScore() const;
	void basicScaleScore(double score);
	
private:
	// The simulation of the actor-based model that owns this variable
	EpochSimulation * lpSimulation;

	// The underlying set of actors
	//const SimulationActorSet * lpActorSet;

	// The current period, in [0, observations - 2]
	int lperiod {};

	// Model parameters (excl B)
	double lbasicScale {};		// Basic scale parameter for the current period
	double lA {};				// Feedback matrix - MatrixXd A(1,1);
	double lG {};				// Diffusion parameters - MatrixXd * G;

	// Bergstrom coefficients (excl B)
	double lAdt {};
	double lQdt {};
	
	// The score for the basic scale parameter for this variable for this period
	double lbasicScaleScore {};

	// The derivative for the basic rate parameter for this variable for
	// this period - only for MLE
	// double lbasicScaleDerivative {};
};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

double SdeSimulation::feedbackParameter() const
{
	return this->lA;
}

double SdeSimulation::feedbackCoefficient() const
{
	return this->lAdt;
}

double SdeSimulation::wienerParameter() const
{
	return this->lG;
}

double SdeSimulation::wienerCoefficient() const 
{
	return this->lQdt;
}

}

#endif /*SDESIMULATION_H_*/

