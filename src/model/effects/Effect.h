/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Effect.h
 *
 * Description: This class contains the declaration of the class Effect.
 *****************************************************************************/

#ifndef EFFECT_H_
#define EFFECT_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class EpochSimulation;
class EffectInfo;
class State;
class Data;
class Cache;


// ----------------------------------------------------------------------------
// Section: Effect class
// ----------------------------------------------------------------------------

/**
 * This class is the base class for all effects. Each effect corresponds to a
 * component in the evaluation or endowment functions, which owns the effect.
 *
 * See the SIENA manual for the mathematical definition of all effects.
 */
class Effect
{
public:
	Effect(const EffectInfo * pEffectInfo);
	virtual ~Effect();

	inline double parameter() const;
	void parameter(double value) ;

	inline int period() const;
	inline Cache * pCache() const;

	const EffectInfo * pEffectInfo() const;

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

private:
	// The effect info object that this effect is created for
	const EffectInfo * lpEffectInfo;

	// The coefficient of this effect in the owner function
	double lparameter {};

	// The period of interest
	int lperiod {};

	Cache * lpCache;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the weight of this effect in the owner function, which is a weighted
 * sum of effects.
 */
double Effect::parameter() const
{
	return this->lparameter;
}


/**
 * Returns the current period.
 */
int Effect::period() const
{
	return this->lperiod;
}


Cache * Effect::pCache() const
{
	return this->lpCache;
}

}

#endif /*EFFECT_H_*/
