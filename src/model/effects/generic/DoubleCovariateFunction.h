/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleCovariateFunction.h
 *
 * Description: This file contains the definition of the
 * DoubleCovariateFunction class.
 *****************************************************************************/

#ifndef DOUBLECOVARIATEFUNCTION_H_
#define DOUBLECOVARIATEFUNCTION_H_

#include "AlterFunction.h"
#include <string>

namespace siena 
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;
class ChangingCovariate;
class BehaviorLongitudinalData;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines an alter function that depends on two covariates or
 * behavior variables.
 */
class DoubleCovariateFunction: public AlterFunction 
{
public:
	DoubleCovariateFunction(std::string covariateName1, std::string covariateName2);
	virtual void initialize(const Data * pData, State * pState, int period,
			Cache * pCache);
	virtual void preprocessEgo(int ego);

protected:
	double firstCovariateValue(int i) const;
	double secondCovariateValue(int i) const;
	int firstCovariateIntValue(int i) const;
	int secondCovariateIntValue(int i) const;
	int firstCovariateN() const;
	int secondCovariateN() const;
	bool firstMissing(int i) const;
	bool secondMissing(int i) const;
	double firstCovariateMinimum() const;
	double secondCovariateMinimum() const;
	double firstCovariateMaximum() const;
	double secondCovariateMaximum() const;

private:
	std::string lFirstCovariateName {};
	std::string lSecondCovariateName {};
	ConstantCovariate * lpFirstConstantCovariate;
	ConstantCovariate * lpSecondConstantCovariate;
	ChangingCovariate * lpFirstChangingCovariate;
	ChangingCovariate * lpSecondChangingCovariate;
	BehaviorLongitudinalData * lpFirstBehaviorData;
    BehaviorLongitudinalData * lpSecondBehaviorData;

	// The current non-centered values of a behavior variable for each actor:
	// These arrays are 0 for covariate-based effects.
	const int * lFirstValues {};
	const int * lSecondValues {};
	int lperiod {};

};

}

#endif /* DOUBLECOVARIATEFUNCTION_H_*/
