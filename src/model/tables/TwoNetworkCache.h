/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TwoNetworkCache.h
 *
 * Description: This file contains the definition of the
 * TwoNetworkCache class.
 *****************************************************************************/


#ifndef TWONETWORKCACHE_H_
#define TWONETWORKCACHE_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class MixedConfigurationTable;
class MixedEgocentricConfigurationTable;
class Network;


// ----------------------------------------------------------------------------
// Section: TwoNetworkCache definition
// ----------------------------------------------------------------------------

/**
 * This class stores varied information regarding a specific ego in a network
 * for repeated use.
 */
class TwoNetworkCache
{
public:
	TwoNetworkCache(const Network * pFirstNetwork,
		const Network *pSecondNetwork);
	virtual ~TwoNetworkCache();

	inline const Network * pFirstNetwork() const;
	inline const Network * pSecondNetwork() const;

	void initialize(int ego);

	inline bool firstOutTieExists(int alter) const;
	inline bool secondOutTieExists(int alter) const;
	inline int firstOutTieValue(int alter) const;
	inline int secondOutTieValue(int alter) const;
	inline int secondOutDegree() const;

	inline MixedEgocentricConfigurationTable * pTwoPathTable() const;
	inline MixedEgocentricConfigurationTable * pInStarTable() const;
	inline MixedEgocentricConfigurationTable * pOutStarTable() const;

private:
	// The first network this cache object is associated with
	const Network * lpFirstNetwork;

	// The second network this cache object is associated with
	const Network * lpSecondNetwork;

	bool loneModeFirstNetwork;
	bool loneModeSecondNetwork;

	// Stores the values of ties from ego in first network
	// to each of the alters who are senders in the second network.
	int * lfirstOutTieValues;

	// Stores the values of ties from ego in first network
	// to each of the alters who are senders in the second network.
	int * lsecondOutTieValues;
	int lsecondOutDegree;

	// The number of two-paths from the ego to each of the alters
	MixedEgocentricConfigurationTable * lpTwoPathTable;

	// The number of in-two-stars from ego and alter to the middle actor.
	MixedEgocentricConfigurationTable * lpInStarTable;

	// The number of out-two-stars from the middle actor to ego and alter.
	MixedEgocentricConfigurationTable * lpOutStarTable;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the first network this cache instance is associated with.
 */
const Network * TwoNetworkCache::pFirstNetwork() const
{
	return this->lpFirstNetwork;
}

/**
 * Returns the second network this cache instance is associated with.
 */
const Network * TwoNetworkCache::pSecondNetwork() const
{
	return this->lpSecondNetwork;
}

/**
 * Indicates if there is a tie from the ego to the given alter in the first
 * network. This method runs in constant time.
 */
bool TwoNetworkCache::firstOutTieExists(int alter) const
{
	return (this->lfirstOutTieValues[alter] >= 1);
}

/**
 * Indicates if there is a tie from the ego to the given alter in the second
 * network. This method runs in constant time.
 */
bool TwoNetworkCache::secondOutTieExists(int alter) const
{
	return (this->lsecondOutTieValues[alter] >= 1);
}


/**
 * Returns the value of the tie from the ego to the given alter in the
 * first network.
 */
int TwoNetworkCache::firstOutTieValue(int alter) const
{
	return this->lfirstOutTieValues[alter];
}


/**
 * Returns the value of the tie from the ego to the given alter in the
 * second network.
 */
int TwoNetworkCache::secondOutTieValue(int alter) const
{
	return this->lsecondOutTieValues[alter];
}

/**
 * Returns the outdegree of the ego in the second network.
 */
int TwoNetworkCache::secondOutDegree() const
{
	return this->lsecondOutDegree;
}

/**
 * Returns the table storing the number of two-paths from the ego to
 * each of the alters.
 */
MixedEgocentricConfigurationTable * TwoNetworkCache::pTwoPathTable() const
{
	return this->lpTwoPathTable;
}

/**
 * Returns the table storing the number of in-two-stars from ego
 * and alter to the middle actor.
 */
MixedEgocentricConfigurationTable * TwoNetworkCache::pInStarTable() const
{
	return this->lpInStarTable;
}

/**
 * Returns the table storing the number of out-two-stars to ego
 * and alter from the middle actor.
 */
MixedEgocentricConfigurationTable * TwoNetworkCache::pOutStarTable() const
{
	return this->lpOutStarTable;
}

}

#endif /* TWONETWORKCACHE_H_ */
