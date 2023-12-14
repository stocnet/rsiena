/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantFunction.h
 *
 * Description: This file contains the definition of the
 * ConstantFunction class.
 *****************************************************************************/

#ifndef CONSTANTFUNCTION_H_
#define CONSTANTFUNCTION_H_

#include <string>
#include "AlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

/**
 * This enumeration defines the possible constants, which may or may not be
 * specific to a dependent variable.
 */
enum ConstantType {VALUE, AVERAGE_IN_DEGREE, AVERAGE_OUT_DEGREE, AVERAGE_RECIP_DEGREE};


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

class ConstantFunction: public AlterFunction
{
public:
	ConstantFunction(double constant);
	ConstantFunction(std::string variableName, ConstantType constantType);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter) const;
	void pFunction(double (* pFunction)(double));

private:
	bool networkConstant() const;

	double lconstant {};

	// The type of the constant that is to be read from the observed
	// data of a dependent variable, or VALUE if a plain constant
	// is to be returned.

	ConstantType lconstantType {};

	// The name of the variable if the constant is to be read from the
	// observed data of a dependent variable.

	std::string lvariableName {};

	// The function to be applied on the constant
	double (* lpFunction)(double);
};

}

#endif /* CONSTANTFUNCTION_H_ */
