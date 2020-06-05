/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: Function.h
 * 
 * Description: This file defines the Function class.
 *****************************************************************************/

#ifndef FUNCTION_H_
#define FUNCTION_H_

#include <vector>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Effect;


// ----------------------------------------------------------------------------
// Section: Function class
// ----------------------------------------------------------------------------

/**
 * This class defines an evaluation or endowment function as a collection of
 * effects.
 */
class Function
{
public:
	virtual ~Function();
	
	void addEffect(Effect * pEffect);
	const std::vector<Effect *> & rEffects() const;

private:
	// A list of effects defining this function
	std::vector<Effect *> leffects;
};

}

#endif /*FUNCTION_H_*/
