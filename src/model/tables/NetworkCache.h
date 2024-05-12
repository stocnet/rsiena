/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkCache.h
 *
 * Description: This file contains the definition of the
 * NetworkCache class.
 *****************************************************************************/


#ifndef NETWORKCACHE_H_
#define NETWORKCACHE_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConfigurationTable;
class EgocentricConfigurationTable;
class Network;


// ----------------------------------------------------------------------------
// Section: NetworkCache definition
// ----------------------------------------------------------------------------

/**
 * This class stores varied information regarding a specific ego in a network
 * for repeated use.
 */
class NetworkCache
{
public:
	NetworkCache(const Network * pNetwork);
	virtual ~NetworkCache();

	inline const Network * pNetwork() const;

	void initialize(int ego);

	inline bool outTieExists(int alter) const;
	inline bool inTieExists(int alter) const;
	inline int outTieValue(int alter) const;
	inline int inTieValue(int alter) const;
	inline int stepTypeValue() const;
	inline void stepTypeSet(int stepType);

	inline EgocentricConfigurationTable * pTwoPathTable() const;
	inline EgocentricConfigurationTable * pReverseTwoPathTable() const;
	inline EgocentricConfigurationTable * pInStarTable() const;
	inline EgocentricConfigurationTable * pOutStarTable() const;
	inline EgocentricConfigurationTable * pCriticalInStarTable() const;
	inline EgocentricConfigurationTable * pRRTable() const;
	inline EgocentricConfigurationTable * pRFTable() const;
	inline EgocentricConfigurationTable * pRBTable() const;
	inline EgocentricConfigurationTable * pFRTable() const;
	inline EgocentricConfigurationTable * pBRTable() const;
	inline ConfigurationTable * pBetweennessTable() const;

private:
	// The network this cache object is associated with
	const Network * lpNetwork;

	bool loneModeNetwork {};

	// Stores the values of ties from ego to each of the alters.
	int * loutTieValues {};

	// Stores the values of ties to ego from each of the alters.
	int * linTieValues {};

	// Stores the values of the stepType (settings model)
	int lstepTypeValue {};

	// The number of two-paths from the ego to each of the alters
	EgocentricConfigurationTable * lpTwoPathTable;

	// The number of two-paths from each of the alters to the ego
	EgocentricConfigurationTable * lpReverseTwoPathTable;

	// The number of in-stars between the ego and each of the alters.
	EgocentricConfigurationTable * lpInStarTable;

	// The number of out-stars between the ego and each of the alters.
	EgocentricConfigurationTable * lpOutStarTable;

	// The number of in-stars <(i,h), (j,h)> between the ego i and each
	// of the alters j, such that there are no two paths i -> h' -> h for
	// h' != j.

	EgocentricConfigurationTable * lpCriticalInStarTable;

	// The number of actors h with reciprocated ties to both i and j.
	EgocentricConfigurationTable * lpRRTable;

	// The number of actors h with a reciprocated tie to i and a tie to j.
	EgocentricConfigurationTable * lpRFTable;

	// The number of actors h with a reciprocated tie to i and a tie from j.
	EgocentricConfigurationTable * lpRBTable;

	// The number of actors h with a tie to i and a reciprocated tie to j.
	EgocentricConfigurationTable * lpFRTable;

	// The number of actors h with a tie from i and a reciprocated tie to j.
	EgocentricConfigurationTable * lpBRTable;

	// The number of non-transitive two-paths through each actor.
	ConfigurationTable * lpBetweennessTable;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network this cache instance is associated with.
 */
const Network * NetworkCache::pNetwork() const
{
	return this->lpNetwork;
}


/**
 * Indicates if there is a tie from the ego to the given alter. This
 * method runs in constant time.
 */
bool NetworkCache::outTieExists(int alter) const
{
	return this->loutTieValues[alter];
}


/**
 * Indicates if there is a tie from the given alter to the ego. This
 * method runs in constant time.
 */
bool NetworkCache::inTieExists(int alter) const
{
	return this->linTieValues[alter];
}


/**
 * Returns the value of the tie from the ego to the given alter.
 */
int NetworkCache::outTieValue(int alter) const
{
	return this->loutTieValues[alter];
}


/**
 * Returns the value of the tie from the given alter to the ego.
 */
int NetworkCache::inTieValue(int alter) const
{
	return this->linTieValues[alter];
}


/**
 * Returns the value of the stepType.
 */
int NetworkCache::stepTypeValue() const
{
	return this->lstepTypeValue;
}

void NetworkCache::stepTypeSet(int stepType)
{
	this->lstepTypeValue = stepType;
}

/**
 * Returns the table storing the number of two-paths from the ego to
 * each of the alters.
 */
EgocentricConfigurationTable * NetworkCache::pTwoPathTable() const
{
	return this->lpTwoPathTable;
}


/**
 * Returns the table storing the number of two-paths from each of the
 * alters to the ego.
 */
EgocentricConfigurationTable * NetworkCache::pReverseTwoPathTable() const
{
	return this->lpReverseTwoPathTable;
}


/**
 * Returns the table storing the number of in-stars between the ego and
 * each of the alters.
 */
EgocentricConfigurationTable * NetworkCache::pInStarTable() const
{
	return this->lpInStarTable;
}


/**
 * Returns the table storing the number of out-stars between the ego and
 * each of the alters.
 */
EgocentricConfigurationTable * NetworkCache::pOutStarTable() const
{
	return this->lpOutStarTable;
}


/**
 * Returns the table storing the number of critical in-stars between the
 * ego and each of the alters. An in-star <(i,h), (j,h)> is critical if
 * there are no two paths i -> h' -> h for h' != j.
 */
EgocentricConfigurationTable * NetworkCache::pCriticalInStarTable() const
{
	return this->lpCriticalInStarTable;
}


/**
 * Returns the table storing the number of actors with reciprocated ties
 * to both i and j.
 */
EgocentricConfigurationTable * NetworkCache::pRRTable() const
{
	return this->lpRRTable;
}


/**
 * Returns the table storing the number of actors with a reciprocated tie
 * to i and a tie to j.
 */
EgocentricConfigurationTable * NetworkCache::pRFTable() const
{
	return this->lpRFTable;
}


/**
 * Returns the table storing the number of actors with a reciprocated tie
 * to i and a tie from j.
 */
EgocentricConfigurationTable * NetworkCache::pRBTable() const
{
	return this->lpRBTable;
}


/**
 * Returns the table storing the number of actors with a tie to i and a
 * reciprocated tie to j.
 */
EgocentricConfigurationTable * NetworkCache::pFRTable() const
{
	return this->lpFRTable;
}


/**
 * Returns the table storing the number of actors with a tie from i and a
 * reciprocated tie to j.
 */
EgocentricConfigurationTable * NetworkCache::pBRTable() const
{
	return this->lpBRTable;
}


/**
 * Returns the table storing the number of non-transitive two-paths
 * through each actor.
 */
ConfigurationTable * NetworkCache::pBetweennessTable() const
{
	return this->lpBetweennessTable;
}

}

#endif /* NETWORKCACHE_H_ */
