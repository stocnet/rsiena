/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Cache.h
 *
 * Description: This file contains the definition of the Cache class.
 *****************************************************************************/

#ifndef CACHE_H_
#define CACHE_H_

#include <map>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class TwoNetworkCache;
class NetworkCache;
class Network;
// ----------------------------------------------------------------------------
// Section: Cache class
// ----------------------------------------------------------------------------

/**
 * This class provides cache objects
 */

class Cache
{
public:
	Cache();
	virtual ~Cache();

	NetworkCache * pNetworkCache(const Network * pNetwork);
	TwoNetworkCache * pTwoNetworkCache(const Network * pFirstNetwork,
		const Network * pSecondNetwork);
	void initialize(int ego);

private:
	std::map<const Network *, NetworkCache *> lnetworkCaches;
	std::map<const Network *, std::map<const Network *, TwoNetworkCache *> >
		ltwoNetworkCaches;
	int lego {};
};

}

#endif /* CACHE_H_ */
