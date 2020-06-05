/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Cache.cpp
 *
 * Description: This file contains the implementation of the class Cache.
 *****************************************************************************/

#include "Cache.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/TwoNetworkCache.h"

using namespace std;

namespace siena
{
// ----------------------------------------------------------------------------
// Section: Constructors and destructors
// ----------------------------------------------------------------------------

/**
 * Creates a cache
 *
 */

Cache::Cache()
{
	this->lego = -1;
}


Cache::~Cache()
{
	clearMap(this->lnetworkCaches, false, true);
	for (map<const Network *, map<const Network *, TwoNetworkCache *> >
			 ::iterator iter = this->ltwoNetworkCaches.begin();
		 iter != this->ltwoNetworkCaches.end();
		 iter++)
	{
		clearMap(iter->second, false, true);
	}
}


NetworkCache * Cache::pNetworkCache(const Network * pNetwork)
{
	NetworkCache * pNetworkCache = 0;
	map<const Network *, NetworkCache *>::iterator iter =
		this->lnetworkCaches.find(pNetwork);

	if (iter != this->lnetworkCaches.end())
	{
		pNetworkCache = iter->second;
	}
	else
	{
		pNetworkCache = new NetworkCache(pNetwork);
		pNetworkCache->initialize(this->lego);
		this->lnetworkCaches[pNetwork] = pNetworkCache;
	}

	return pNetworkCache;
}

TwoNetworkCache * Cache::pTwoNetworkCache(const Network * pFirstNetwork,
	const Network * pSecondNetwork)
{
	TwoNetworkCache * pTwoNetworkCache = 0;
	map<const Network *,
		map <const Network *, TwoNetworkCache *> >::iterator iter =
		this->ltwoNetworkCaches.find(pFirstNetwork);

	if (iter != this->ltwoNetworkCaches.end())
	{
		map <const Network *, TwoNetworkCache *> cacheMap = iter->second;
		map <const Network *, TwoNetworkCache *>::iterator iter2 =
			cacheMap.find(pSecondNetwork);
		if (iter2 != cacheMap.end())
		{
			pTwoNetworkCache = iter2->second;
		}
	}
	if (!pTwoNetworkCache)
	{
		pTwoNetworkCache = new TwoNetworkCache(pFirstNetwork, pSecondNetwork);
		pTwoNetworkCache->initialize(this->lego);

		this->ltwoNetworkCaches[pFirstNetwork][pSecondNetwork] =
			pTwoNetworkCache;
	}

	return pTwoNetworkCache;
}

void Cache::initialize(int ego)
{
	this->lego = ego;

	for (map<const Network *, NetworkCache *>::iterator iter =
			this->lnetworkCaches.begin();
		iter != this->lnetworkCaches.end();
		iter++)
	{
		NetworkCache * pNetworkCache = iter->second;
		pNetworkCache->initialize(ego);
	}

	for (map<const Network *, map<const Network *, TwoNetworkCache *> >
			 ::iterator iter = this->ltwoNetworkCaches.begin();
		 iter != this->ltwoNetworkCaches.end();
		 iter++)
	{
		map<const Network *, TwoNetworkCache *> cacheMap = iter->second;
		for (map<const Network *, TwoNetworkCache *>::iterator iter2 =
				 cacheMap.begin();
			 iter2 != cacheMap.end();
			 iter2++)
		{
			TwoNetworkCache * pTwoNetworkCache = iter2->second;
			pTwoNetworkCache->initialize(ego);
		}
	}
}

}
