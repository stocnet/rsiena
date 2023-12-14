/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkEffect.h
 *
 * Description: This file contains the definition of the
 * NetworkEffect class.
 *****************************************************************************/

#ifndef NETWORKEFFECT_H_
#define NETWORKEFFECT_H_

#include "Effect.h"
#include <utility>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class NetworkLongitudinalData;
class ConfigurationTable;
class NetworkCache;
class Cache;


// ----------------------------------------------------------------------------
// Section: NetworkEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all network effects.
 */
class NetworkEffect : public Effect
{
	friend class NetworkInteractionEffect;
	friend class EffectFactory;
	friend class BothDegreesEffect;

public:
	NetworkEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void initialize(const Data * pData, State * pState,
			State * pSimulatedState, int period, Cache * pCache);

	inline const Network * pNetwork() const;
	inline const NetworkLongitudinalData * pData() const;

	virtual void preprocessEgo(int ego);
	inline int ego() const;

	/**
	 * Assuming that the ego would flip the tie to the given actor,
	 * this method calculates the change in the statistic corresponding
	 * to this effect. The method has to be overriden by all concrete
	 * effect classes.
	 */
	virtual double calculateContribution(int alter) const = 0;

	virtual double evaluationStatistic();
	virtual std::pair<double, double * > evaluationStatistic(bool needActorStatistics);
	virtual double endowmentStatistic(Network * pLostTieNetwork);
	virtual std::pair<double, double * > endowmentStatistic(Network * pLostTieNetwork, bool needActorStatistics);
	virtual double creationStatistic(Network * pGainedTieNetwork);
	virtual std::pair<double, double * > creationStatistic(Network * pGainedTieNetwork, bool needActorStatistics);

	virtual bool egoEffect() const;

protected:
	int n() const;
	virtual double statistic(const Network * pSummationTieNetwork);
	virtual std::pair<double, double * > statistic(const Network * pSummationTieNetwork, bool needActorStatistics);
	virtual void initializeStatisticCalculation();
	virtual void onNextEgo(int ego);
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);
	virtual double tieStatistic(int alter);
	virtual void cleanupStatisticCalculation();
	inline int settingsStepType() const;

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
	bool inTieExists(int alter) const;
	bool outTieExists(int alter) const;

private:
	// The network this effect is associated with
	const Network * lpNetwork;

	// The observed network data underlying the network variable.
	const NetworkLongitudinalData * lpNetworkData;

	NetworkCache * lpNetworkCache;
	int lego {};

	// Stores the values of the stepType (settings model)
	int lstepTypeVal {};

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
 * Returns the network this effect is associated with.
 */
const Network * NetworkEffect::pNetwork() const
{
	return this->lpNetwork;
}


/**
 * Returns the observed network data corresponding to the network variable of
 * this effect.
 */
const NetworkLongitudinalData * NetworkEffect::pData() const
{
	return this->lpNetworkData;
}


/**
 * Returns the current ego, whose tie flip contributions are to be calculated.
 */
int NetworkEffect::ego() const
{
	return this->lego;
}


/**
 * Returns the stepType of the current ministep.
 */
int NetworkEffect::settingsStepType() const
{
	return this->lstepTypeVal;
}


/**
 * Returns the table storing the number of two-paths from the ego to
 * each of the alters.
 */
inline ConfigurationTable * NetworkEffect::pTwoPathTable() const
{
	return this->lpTwoPathTable;
}


/**
 * Returns the table storing the number of two-paths from each of the
 * alters to the ego.
 */
inline ConfigurationTable * NetworkEffect::pReverseTwoPathTable() const
{
	return this->lpReverseTwoPathTable;
}


/**
 * Returns the table storing the number of in-stars between the ego and
 * each of the alters.
 */
inline ConfigurationTable * NetworkEffect::pInStarTable() const
{
	return this->lpInStarTable;
}


/**
 * Returns the table storing the number of out-stars between the ego and
 * each of the alters.
 */
inline ConfigurationTable * NetworkEffect::pOutStarTable() const
{
	return this->lpOutStarTable;
}


/**
 * Returns the table storing the number of critical in-stars between the
 * ego and each of the alters. An in-star <(i,h), (j,h)> is critical if
 * there are no two paths i -> h' -> h for h' != j.
 */
inline ConfigurationTable * NetworkEffect::pCriticalInStarTable() const
{
	return this->lpCriticalInStarTable;
}


/**
 * Returns the table storing the number of actors with reciprocated ties
 * to both i and j.
 */
inline ConfigurationTable * NetworkEffect::pRRTable() const
{
	return this->lpRRTable;
}


/**
 * Returns the table storing the number of actors with a reciprocated tie
 * to i and a tie to j.
 */
inline ConfigurationTable * NetworkEffect::pRFTable() const
{
	return this->lpRFTable;
}


/**
 * Returns the table storing the number of actors with a reciprocated tie
 * to i and a tie from j.
 */
inline ConfigurationTable * NetworkEffect::pRBTable() const
{
	return this->lpRBTable;
}


/**
 * Returns the table storing the number of actors with a tie to i and a
 * reciprocated tie to j.
 */
inline ConfigurationTable * NetworkEffect::pFRTable() const
{
	return this->lpFRTable;
}


/**
 * Returns the table storing the number of actors with a tie from i and a
 * reciprocated tie to j.
 */
inline ConfigurationTable * NetworkEffect::pBRTable() const
{
	return this->lpBRTable;
}

}

#endif /*NETWORKEFFECT_H_*/
