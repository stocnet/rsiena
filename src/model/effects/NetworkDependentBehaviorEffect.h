/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * NetworkDependentBehaviorEffect class.
 *****************************************************************************/

#ifndef NETWORKDEPENDENTBEHAVIOREFFECT_H_
#define NETWORKDEPENDENTBEHAVIOREFFECT_H_

#include "BehaviorEffect.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class ConfigurationTable;
class NetworkCache;

// ----------------------------------------------------------------------------
// Section: NetworkDependentBehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects depending on some network variable.
 */
class NetworkDependentBehaviorEffect : public BehaviorEffect
{
public:
	NetworkDependentBehaviorEffect(const EffectInfo * pEffectInfo);
	// for gmom:
	NetworkDependentBehaviorEffect(const EffectInfo * pEffectInfo, const bool simulatedState);
	virtual ~NetworkDependentBehaviorEffect();

	virtual void initialize(const Data * pData,
			State * pState, int period, Cache * pCache);
	virtual void initialize(const Data *pData,
			State *pState, State *pSimulatedState, int period, Cache *pCache);

protected:
	inline const Network * pNetwork() const;
	double totalAlterValue(int i) const;
	double totalInAlterValue(int i) const;
	int numberAlterHigher(int i) const;
	int numberAlterLower(int i) const;
	int numberAlterEqual(int i) const;
	int numberAlterHigherPop(int i) const;
	int numberAlterLowerPop(int i) const;
	int numberAlterEqualPop(int i) const;
	virtual void preprocessEgo(int ego);

	inline ConfigurationTable * pTwoPathTable() const;
	inline ConfigurationTable * pReverseTwoPathTable() const;
	inline ConfigurationTable * pInStarTable() const;
	inline ConfigurationTable * pOutStarTable() const;
	inline ConfigurationTable * pCriticalInStarTable() const;
	inline ConfigurationTable * pRRTable() const;
	inline ConfigurationTable * pRFTable() const;
	inline ConfigurationTable * pRBTable() const;
	inline ConfigurationTable * pFRTable() const;
	inline ConfigurationTable * pBRTable() const;

private:
	//! If `1` value(), missing() and similarity() returns the simulated value
	//! (if the covariate is a behavior) or the observed value at the end of the
	//! period.
	const int lSimulatedOffset {};

	// The network this effect is interacting with
	const Network * lpNetwork;
	// total out- and in-alter values
	double * ltotalAlterValues {};
	double * ltotalInAlterValues {};
	// number of higher, lower, and equal alter values
	int * lnumberAlterHigher {};
	int * lnumberAlterLower {};
	int * lnumberAlterEqual {};
	// and weighted by alter indegrees
	int * lnumberAlterHigherPop {};
	int * lnumberAlterLowerPop {};
	int * lnumberAlterEqualPop {};

	NetworkCache * lpNetworkCache;

	// The number of two-paths from the ego to each of the alters
	ConfigurationTable * lpTwoPathTable;

	// The number of two-paths from each of the alters to the ego
	ConfigurationTable * lpReverseTwoPathTable;

	// The number of in-stars between the ego and each of the alters.
	ConfigurationTable * lpInStarTable;

	// The number of out-stars between the ego and each of the alters.
	ConfigurationTable * lpOutStarTable;

	// The number of in-stars <(i,h), (j,h)> between the ego i and each
	// of the alters j, such that there are no two paths i -> h' -> h for
	// h' != j.

	ConfigurationTable * lpCriticalInStarTable;

	// The number of actors h with reciprocated ties to both i and j.
	ConfigurationTable * lpRRTable;

	// The number of actors h with a reciprocated tie to i and a tie to j.
	ConfigurationTable * lpRFTable;

	// The number of actors h with a reciprocated tie to i and a tie from j.
	ConfigurationTable * lpRBTable;

	// The number of actors h with a tie to i and a reciprocated tie to j.
	ConfigurationTable * lpFRTable;

	// The number of actors h with a tie from i and a reciprocated tie to j.
	ConfigurationTable * lpBRTable;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network this effect is interacting with.
 */
const Network * NetworkDependentBehaviorEffect::pNetwork()
	const
{
	return this->lpNetwork;
}

/**
 * Returns the table storing the number of two-paths from the ego to
 * each of the alters.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pTwoPathTable() const
{
	return this->lpTwoPathTable;
}


/**
 * Returns the table storing the number of two-paths from each of the
 * alters to the ego.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pReverseTwoPathTable() const
{
	return this->lpReverseTwoPathTable;
}


/**
 * Returns the table storing the number of in-stars between the ego and
 * each of the alters.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pInStarTable() const
{
	return this->lpInStarTable;
}


/**
 * Returns the table storing the number of out-stars between the ego and
 * each of the alters.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pOutStarTable() const
{
	return this->lpOutStarTable;
}


/**
 * Returns the table storing the number of critical in-stars between the
 * ego and each of the alters. An in-star <(i,h), (j,h)> is critical if
 * there are no two paths i -> h' -> h for h' != j.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pCriticalInStarTable() const
{
	return this->lpCriticalInStarTable;
}


/**
 * Returns the table storing the number of actors with reciprocated ties
 * to both i and j.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pRRTable() const
{
	return this->lpRRTable;
}


/**
 * Returns the table storing the number of actors with a reciprocated tie
 * to i and a tie to j.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pRFTable() const
{
	return this->lpRFTable;
}


/**
 * Returns the table storing the number of actors with a reciprocated tie
 * to i and a tie from j.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pRBTable() const
{
	return this->lpRBTable;
}


/**
 * Returns the table storing the number of actors with a tie to i and a
 * reciprocated tie to j.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pFRTable() const
{
	return this->lpFRTable;
}


/**
 * Returns the table storing the number of actors with a tie from i and a
 * reciprocated tie to j.
 */
inline ConfigurationTable * NetworkDependentBehaviorEffect::pBRTable() const
{
	return this->lpBRTable;
}

}

#endif /*NETWORKDEPENDENTBEHAVIOREFFECT_H_*/
