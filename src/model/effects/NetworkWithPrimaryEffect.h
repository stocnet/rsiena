/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkWithPrimaryEffect.h
 *
 * Description: This file contains the definition of the
 * NetworkWithPrimaryEffect class.
 *****************************************************************************/

#ifndef NETWORKWITHPRIMARYEFFECT_H_
#define NETWORKWITHPRIMARYEFFECT_H_

#include <string>
#include "NetworkEffect.h"
#include "utils/NamedObject.h"
using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class Cache;

// ----------------------------------------------------------------------------
// Section: NetworkWithPrimaryEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all network effects depending on the primary setting.
 */
class NetworkWithPrimaryEffect : public NetworkEffect
{
public:
	NetworkWithPrimaryEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual ~NetworkWithPrimaryEffect();
	virtual void preprocessEgo(int ego);

protected:
// calculate primary setting properties:
	void primaryProperties(int ego, const Network * pNetwork);
// degree of ego in primary setting, after preprocessEgo(int ego):
	int primaryDegree() const;
	bool inPrimarySetting(int alter) const;

private:	
	// vector of indicators of primary setting:
	bool * lprimary {};
	// egos primary setting degree:
	int lprimDegree {};
};
}

#endif /*NETWORKWITHPRIMARYEFFECT_H_*/
