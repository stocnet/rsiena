/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EffectValueTable.h
 *
 * Description: This file contains the definition of the
 * EffectValueTable class.
 *****************************************************************************/

#ifndef EFFECTVALUETABLE_H_
#define EFFECTVALUETABLE_H_

namespace siena
{

/**
 * This class provides a look-up table supporting effective calculation of
 * structural effects on rate functions. The value of the effect for an
 * argument <i>i</i> in [0, <i>n</i>) is defined as
 * exp(alpha <i>f</i>(<i>i</i>)), where alpha is a parameter associated
 * with this effect and <i>f</i> is an arbitrary function. Since calculating
 * exponentials is expensive, the values are stored for each <i>i</i> for
 * later reuse as long as the parameter remains unchanged.
 */
class EffectValueTable
{
public:
	EffectValueTable(int n, double (* pFunction)(int));
	virtual ~EffectValueTable();

	double parameter() const;
	void parameter(double value);
	double value(int i);

private:
	// The function to be applied
	double (* lpFunction)(int);

	// The look-up table for effect values
	double * lvalues {};

	// Here we remember the parameter value whenever we store a value
	// in the table.

	double * lparameterValues {};

	// The actual value of the effect parameter
	double lparameter {};
};

}

#endif /*EFFECTVALUETABLE_H_*/
