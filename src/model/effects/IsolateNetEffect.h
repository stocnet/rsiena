/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IsolateNetEffect.h
 *
 * Description: This file contains the definition of the
 * IsolateNetEffect class.
 *****************************************************************************/

#ifndef ISOLATENETEFFECT_H_
#define ISOLATENETEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * This class defines the network-isolate effect defined by
 * s_i(x)= I{x_{i+}=x_{+i}=0} 
 * The corresponding statistic is
 * the number of isolate.
 */
class IsolateNetEffect : public NetworkEffect
{
public:
	IsolateNetEffect(const EffectInfo * pEffectInfo);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);

};

}

#endif /*ISOLATENETEFFECT_H_*/
