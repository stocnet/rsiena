/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EffectFactory.h
 *
 * Description: This file contains the definition of the
 * EffectFactory class.
 *****************************************************************************/

#ifndef EFFECTFACTORY_H_
#define EFFECTFACTORY_H_

#include <map>
#include <string>
#include <model/EffectInfo.h>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Data;
class Effect;
class EffectInfo;


// ----------------------------------------------------------------------------
// Section: EffectFactory class
// ----------------------------------------------------------------------------

/**
 * Produces concrete effects of the Effect class hierarchy from generic
 * effect descriptor objects of class EffectInfo.
 */
class EffectFactory
{
public:
	static const std::string gmmGroup(const EffectInfo* pEffectInfo);

	EffectFactory(const Data * pData);

	Effect * createEffect(const EffectInfo * pEffectInfo) const;

private:
	const Data * lpData;
	static const std::map<const std::string, const std::string> GMM_GROUPS;
	static std::map<const std::string, const std::string> init_groups();
};

}

#endif /*EFFECTFACTORY_H_*/
