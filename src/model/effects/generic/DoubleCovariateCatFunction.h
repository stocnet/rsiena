/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleCovariateCatFunction.h
 *
 * Description: This file contains the definition of the
 * DoubleCovariateCatFunction class.
 *****************************************************************************/

#ifndef DOUBLECOVARIATECATFUNCTION_H_
#define DOUBLECOVARIATECATFUNCTION_H_

#include "DoubleCovariateFunction.h"
#include <string>


namespace siena 
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

//class ConstantCovariate;
//class ChangingCovariate;
//class BehaviorLongitudinalData;
class Network;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines an alter function that depends on two covariates or
 * behavior variables, and used specifically for sameXVInPop. 
 * Includes a bit of NetworkAlterFunction.h;
 * if more of NetworkAlterFunction is needed, that can be extended.
 */
class DoubleCovariateCatFunction: public DoubleCovariateFunction
{
public:
	DoubleCovariateCatFunction(std::string covariateName1, std::string covariateName2,
						std::string networkName, double parameter, bool excludeMissing,
						bool byTies);
	virtual ~DoubleCovariateCatFunction();
	virtual void initialize(const Data * pData, State * pState, int period,
			Cache * pCache);
	virtual void preprocessEgo(int ego);	
	int numberCovariateTies(int a, int b) const;
	virtual double value(int alter) const;

protected:
	inline const Network * pNetwork() const;
	int firstCovariateNumbers(int a) const;
	int secondCovariateNumbers(int b) const;
	int numberCovariateTies(int b) const;

private:
	const Network * lpNetwork;
	std::string lnetworkName {};

//	std::string lFirstCovariateName {};
//	std::string lSecondCovariateName {};
//	ConstantCovariate * lpFirstConstantCovariate;
//	ConstantCovariate * lpSecondConstantCovariate;
//	ChangingCovariate * lpFirstChangingCovariate;
//	ChangingCovariate * lpSecondChangingCovariate;
//	BehaviorLongitudinalData * lpFirstBehaviorData;
//   BehaviorLongitudinalData * lpSecondBehaviorData;
	bool lexcludeMissing {};
	bool lbyTies {};

	// An array of counts of tie values
	// between nodes with ego's covariate value and alters with specific covariate values
	int * lpNumberTieValues {};
	int * lpTotalCovariateCombinations {};
	int * lpFirstCovariateNumbers {};
	int * lpSecondCovariateNumbers {};
	bool lroot {}; // should the square root be taken?
	bool laverage {}; // should the average be used?
	int lSecondMax {};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

const Network * DoubleCovariateCatFunction::pNetwork() const
{
	return this->lpNetwork;
}

}

#endif /* DOUBLECOVARIATECATFUNCTION_H_ */
