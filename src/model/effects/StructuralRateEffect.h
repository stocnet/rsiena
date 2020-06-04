/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: StructuralRateEffect.h
 *
 * Description: This file contains the definition of the
 * StructuralRateEffect class.
 *****************************************************************************/

#ifndef STRUCTURALRATEEFFECT_H_
#define STRUCTURALRATEEFFECT_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

/**
 * Available types of structural rate effects.
 */
enum StructuralRateEffectType
{
	OUT_DEGREE_RATE,
	IN_DEGREE_RATE,
	RECIPROCAL_DEGREE_RATE,
	INVERSE_OUT_DEGREE_RATE,
	LOG_OUT_DEGREE_RATE
};


// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class NetworkVariable;
class EffectValueTable;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Encapsulates the information necessary for calculating the contributions
 * of a structural rate effect. This includes the effect type (out-degree,
 * in-degree, reciprocal degree, and inverse out-degree), the network variable
 * the effect depends on (i.e. the variable used to access the degrees
 * of actors), and the statistical parameter of the effect.
 */
class StructuralRateEffect
{
public:
	StructuralRateEffect(const NetworkVariable * pVariable,
		StructuralRateEffectType type, double parameter);
	virtual ~StructuralRateEffect();

	double value(int i) const;
	void parameter(double parameterValue);
	double parameter() const;

private:
	// The network variable this effect depends on
	const NetworkVariable * lpVariable;

	// The type of the effect
	StructuralRateEffectType ltype;

	// A table for efficient calculation of contributions. If two actors have
	// the same degrees, then this effect contributes equally for those actors.
	// The table ensures that we don't calculate the same contribution twice.

	EffectValueTable * lpTable;
};

}

#endif /* STRUCTURALRATEEFFECT_H_ */
