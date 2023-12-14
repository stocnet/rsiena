/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DiffusionEffectValueTable.h
 *
 * Description: This file contains the definition of the
 * DiffusionEffectValueTable class.
 *****************************************************************************/

#ifndef DIFFUSIONEFFECTVALUETABLE_H_
#define DIFFUSIONEFFECTVALUETABLE_H_

namespace siena
{

/**
 * This class provides a look-up table supporting effective calculation of
 * diffusion effects on rate functions. The value of the effect for an
 * argument <i>i</i> in [0, <i>n</i>) is defined as
 * exp(alpha <i>f</i>(<i>i</i>)), where alpha is a parameter associated
 * with this effect and <i>f</i> is an arbitrary function. Since calculating
 * exponentials is expensive, the values are stored for each <i>i</i> for
 * later reuse as long as the parameter remains unchanged.
 */
class DiffusionEffectValueTable
{
public:
	DiffusionEffectValueTable(int numeratorRange, int denominatorRange);
	virtual ~DiffusionEffectValueTable();

	double parameter() const;
	void parameter(double value);
	double value(int i);
	double value(int numerator, int denominator);

private:
	// The look-up table for effect values
	double * lvalues {};

	// Here we remember the parameter value whenever we store a value
	// in the table.

	double * lparameterValues {};

	// The actual value of the effect parameter
	double lparameter {};

	// The range of the denominator;
	int ldenominatorRange {};

	// The range of the numerator;
	int lnumeratorRange {};

};

}

#endif /*DIFFUSIONEFFECTVALUETABLE_H_*/
