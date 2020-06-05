/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Data.cpp
 *
 * Description: This file contains the implementation of the Data class.
 *****************************************************************************/

#include "Data.h"
#include "utils/Utils.h"
#include "utils/NamedObject.h"
#include "data/ActorSet.h"
#include "data/NetworkLongitudinalData.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ContinuousLongitudinalData.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/ExogenousEvent.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructor area
// ----------------------------------------------------------------------------

/**
 * Constructs an empty data object for the given number of observations.
 */
Data::Data(int observationCount)
{
	this->lobservationCount = observationCount;

	// Create empty sets of exogenous events.
//	for (int i = 0; i < observationCount - 1; i++)

	for (int i = 0; i < observationCount; i++)
	{
		this->lexogenousEvents.push_back(new EventSet());
	}
}


/**
 * Deallocates the data.
 */
Data::~Data()
{
	// Activity indicators

	for (unsigned i = 0; i < this->lactorSets.size(); i++)
	{
		const ActorSet * pActorSet = this->lactorSets[i];
		bool ** active = this->lactive[pActorSet];

		for (int i = 0; i < pActorSet->n(); i++)
		{
			delete[] active[i];
		}

		delete[] active;
	}

	this->lactive.clear();

	// Various vectors

	deallocateVector(this->ldependentVariableData);
	deallocateVector(this->lconstantCovariates);
	deallocateVector(this->lchangingCovariates);
	deallocateVector(this->lconstantDyadicCovariates);
	deallocateVector(this->lchangingDyadicCovariates);
	deallocateVector(this->lactorSets);
	deallocateVector(this->lnetworkConstraints);

	// Delete the exogenous events.

//	for (int i = 0; i < this->lobservationCount - 1; i++)
	for (int i = 0; i < this->lobservationCount; i++)
	{
		while (!this->lexogenousEvents[i]->empty())
		{
			ExogenousEvent * pEvent = *this->lexogenousEvents[i]->begin();
			this->lexogenousEvents[i]->erase(
				this->lexogenousEvents[i]->begin());
			delete pEvent;
		}

		delete this->lexogenousEvents[i];
	}

	this->lexogenousEvents.clear();
}


// ----------------------------------------------------------------------------
// Section: Various data object factory methods and accessors
// ----------------------------------------------------------------------------

/**
 * Creates a new data object for storing observations of a two-mode network.
 * @param[in] name the name of the network variable
 * @param[in] pSenders the set of actors acting as senders of ties
 * @param[in] pReceivers the set of actors acting as receivers of ties
 */
NetworkLongitudinalData * Data::createNetworkData(std::string name,
	const ActorSet * pSenders,
	const ActorSet * pReceivers)
{
	NetworkLongitudinalData * pNetworkData =
		new NetworkLongitudinalData(this->ldependentVariableData.size(),
			name,
			pSenders,
			pReceivers,
			this->lobservationCount,
			false);
	this->ldependentVariableData.push_back(pNetworkData);
	return pNetworkData;
}


/**
 * Creates a new data object for storing observations of a one-mode network.
 * @param[in] name the name of the network variable
 * @param[in] pActors the set of actors of the network
 */
OneModeNetworkLongitudinalData * Data::createOneModeNetworkData(
	std::string name, const ActorSet * pActors)
{
	OneModeNetworkLongitudinalData * pNetworkData =
		new OneModeNetworkLongitudinalData(this->ldependentVariableData.size(),
			name,
			pActors,
			this->lobservationCount);
	this->ldependentVariableData.push_back(pNetworkData);
	return pNetworkData;
}

OneModeNetworkLongitudinalData * Data::createOneModeSimNetworkData(
	std::string name, const ActorSet * pActors)
{
	OneModeNetworkLongitudinalData * pNetworkData =
		new OneModeNetworkLongitudinalData(this->lsimVariableData.size(),
			name, pActors, this->lobservationCount);
	this->lsimVariableData.push_back(pNetworkData);
	return pNetworkData;
}


/**
 * Creates a new data object for storing observations of a behavior variable
 * for the given set of actors.
 * @param[in] name the name of the behavior variable
 */
BehaviorLongitudinalData * Data::createBehaviorData(std::string name,
	const ActorSet * pActorSet)
{
	// The index of the longidutinal data object
	BehaviorLongitudinalData * pBehaviorData =
		new BehaviorLongitudinalData(this->ldependentVariableData.size(),
			name,
			pActorSet,
			this->lobservationCount);
	this->ldependentVariableData.push_back(pBehaviorData);
	return pBehaviorData;
}


/**
 * Creates a new data object for storing observations of a continuous behavior
 * variable for the given set of actors.
 * @param[in] name the name of the continuous behavior variable
 */
ContinuousLongitudinalData * Data::createContinuousData(std::string name,
	const ActorSet * pActorSet)
{
	// The index of the longitudinal data object
	ContinuousLongitudinalData * pContinuousData =
		new ContinuousLongitudinalData(this->ldependentVariableData.size(),
			name,
			pActorSet,
			this->lobservationCount);
	this->ldependentVariableData.push_back(pContinuousData);
	return pContinuousData;
}


/**
 * Creates a new data object for storing the values of a constant covariate
 * for the given set of actors.
 * @param[in] name the name of the covariate
 */
ConstantCovariate * Data::createConstantCovariate(std::string name,
	const ActorSet * pActorSet)
{
	ConstantCovariate * pCovariate = new ConstantCovariate(name, pActorSet);
	this->lconstantCovariates.push_back(pCovariate);
	return pCovariate;
}


/**
 * Creates a new data object for storing the values of a changing covariate
 * for the given set of actors.
 * @param[in] name the name of the covariate
 */
ChangingCovariate * Data::createChangingCovariate(std::string name,
	const ActorSet * pActorSet)
{
	ChangingCovariate * pCovariate =
		new ChangingCovariate(name, pActorSet, this->lobservationCount-1);
	this->lchangingCovariates.push_back(pCovariate);
	return pCovariate;
}


/**
 * Creates a new data object for storing the values of a constant dyadic
 * covariate for the given pair of actor sets.
 * @param[in] name the name of the covariate
 */
ConstantDyadicCovariate * Data::createConstantDyadicCovariate(std::string name,
	const ActorSet * pFirstActorSet,
	const ActorSet * pSecondActorSet)
{
	ConstantDyadicCovariate * pCovariate =
		new ConstantDyadicCovariate(name, pFirstActorSet, pSecondActorSet);
	this->lconstantDyadicCovariates.push_back(pCovariate);
	return pCovariate;
}


/**
 * Creates a new data object for storing the values of a changing dyadic
 * covariate for the given pair of actor sets.
 * @param[in] name the name of the covariate
 */
ChangingDyadicCovariate * Data::createChangingDyadicCovariate(std::string name,
	const ActorSet * pFirstActorSet,
	const ActorSet * pSecondActorSet)
{
	ChangingDyadicCovariate * pCovariate =
		new ChangingDyadicCovariate(name,
			pFirstActorSet,
			pSecondActorSet,
			this->lobservationCount);
	this->lchangingDyadicCovariates.push_back(pCovariate);
	return pCovariate;
}


/**
 * Returns the collection of data objects for dependent variables.
 */
const std::vector<LongitudinalData *> & Data::rDependentVariableData() const
{
	return this->ldependentVariableData;
}

const std::vector<LongitudinalData *> & Data::rSimVariableData() const
{
	return this->lsimVariableData;
}

/**
 * Returns the collection of constant covariates.
 */
const std::vector<ConstantCovariate *> & Data::rConstantCovariates() const
{
	return this->lconstantCovariates;
}


/**
 * Returns the collection of changing covariates.
 */
const std::vector<ChangingCovariate *> & Data::rChangingCovariates() const
{
	return this->lchangingCovariates;
}


/**
 * Returns the collection of constant dyadic covariates.
 */
const std::vector<ConstantDyadicCovariate *> &
	Data::rConstantDyadicCovariates() const
{
	return this->lconstantDyadicCovariates;
}


/**
 * Returns the collection of changing dyadic covariates.
 */
const std::vector<ChangingDyadicCovariate *> &
	Data::rChangingDyadicCovariates() const
{
	return this->lchangingDyadicCovariates;
}


// ----------------------------------------------------------------------------
// Section: Actor sets
// ----------------------------------------------------------------------------

/**
 * Creates a new set of actors of the given size and returns the resulting
 * set.
 */
const ActorSet * Data::createActorSet(std::string name, int n)
{
	ActorSet * pActorSet = new ActorSet(name, n);

	this->lactorSets.push_back(pActorSet);

	// Allocate and initialize activity indicators.

	this->lactive[pActorSet] = new bool * [n];

	for (int i = 0; i < n; i++)
	{
		this->lactive[pActorSet][i] = new bool[this->lobservationCount];

		// All actors are active by default

		for (int k = 0; k < this->lobservationCount; k++)
		{
			this->lactive[pActorSet][i][k] = true;
		}
	}

	return pActorSet;
}


/**
 * Returns the collection of actor sets.
 */
const std::vector<const ActorSet *> & Data::rActorSets() const
{
	return this->lactorSets;
}


// ----------------------------------------------------------------------------
// Section: Lookup by name
// ----------------------------------------------------------------------------

/**
 * Searches the given vector for an object with the given name.
 * @return an object with the given name or 0 if nothing was found
 */
template<class T>
T * findNamedObject(std::string name,
	const std::vector<T *> & rVector)
{
	T * pObject = 0;

	for (unsigned i = 0; i < rVector.size() && !pObject; i++)
	{
		if (rVector[i]->name() == name)
		{
			pObject = rVector[i];
		}
	}

	return pObject;
}


/**
 * Returns the actor set with the given name.
 */
const ActorSet * Data::pActorSet(std::string name) const
{
	return findNamedObject(name, this->lactorSets);
}


/**
 * Returns the longitudinal network data with the given name.
 */
NetworkLongitudinalData * Data::pNetworkData(std::string name) const
{
	return dynamic_cast<NetworkLongitudinalData *>(
		findNamedObject(name, this->ldependentVariableData));
}

NetworkLongitudinalData * Data::pSimNetworkData(std::string name) const
{
	return dynamic_cast<NetworkLongitudinalData *>(
		findNamedObject(name, this->lsimVariableData));
}

/**
 * Returns the longitudinal one-mode network data with the given name.
 */
OneModeNetworkLongitudinalData * Data::pOneModeNetworkData(std::string name) const
{
	return dynamic_cast<OneModeNetworkLongitudinalData *>(
		findNamedObject(name, this->ldependentVariableData));
}

OneModeNetworkLongitudinalData * Data::pOneModeSimNetworkData(std::string name) const
{
	return dynamic_cast<OneModeNetworkLongitudinalData *>(
		findNamedObject(name, this->lsimVariableData));
}


/**
 * Returns the longitudinal behavior data with the given name.
 */
BehaviorLongitudinalData * Data::pBehaviorData(std::string name) const
{
	return dynamic_cast<BehaviorLongitudinalData *>(
		findNamedObject(name, this->ldependentVariableData));
}


/**
 * Returns the longitudinal continuous behavior data with the given name.
 */
ContinuousLongitudinalData * Data::pContinuousData(std::string name) const
{
	return dynamic_cast<ContinuousLongitudinalData *>(
		findNamedObject(name, this->ldependentVariableData));
}


/**
 * Returns the constant covariate with the given name.
 */
ConstantCovariate * Data::pConstantCovariate(std::string name) const
{
	return findNamedObject(name, this->lconstantCovariates);
}


/**
 * Returns the changing covariate with the given name.
 */
ChangingCovariate * Data::pChangingCovariate(std::string name) const
{
	return findNamedObject(name, this->lchangingCovariates);
}


/**
 * Returns the constant dyadic covariate with the given name.
 */
ConstantDyadicCovariate * Data::pConstantDyadicCovariate(std::string name)
	const
{
	return findNamedObject(name, this->lconstantDyadicCovariates);
}


/**
 * Returns the changing dyadic covariate with the given name.
 */
ChangingDyadicCovariate * Data::pChangingDyadicCovariate(std::string name)
	const
{
	return findNamedObject(name, this->lchangingDyadicCovariates);
}


// ----------------------------------------------------------------------------
// Section: Activity and composition change.
// ----------------------------------------------------------------------------

/**
 * Sets if the given actor of the given set is active at the given observation.
 */
void Data::active(const ActorSet * pActorSet,
	int actor,
	int observation,
	bool flag)
{
	this->lactive[pActorSet][actor][observation] = flag;
}


/**
 * Returns if the given actor of the given set is active at the given
 * observation.
 */
bool Data::active(const ActorSet * pActorSet, int actor, int observation)
{
	return this->lactive[pActorSet][actor][observation];
}


/**
 * Adds an exogenous composition change event, specifying that the given
 * actor of the given set joins at the given time of the given period.
 */
void Data::addJoiningEvent(int period,
	const ActorSet * pActorSet,
	int actor,
	double time)
{
	this->lexogenousEvents[period]->insert(
		new ExogenousEvent(pActorSet, actor, time, JOINING));
}


/**
 * Adds an exogenous composition change event, specifying that the given
 * actor of the given set leaves at the given time of the given period.
 */
void Data::addLeavingEvent(int period,
	const ActorSet * pActorSet,
	int actor,
	double time)
{
	this->lexogenousEvents[period]->insert(
		new ExogenousEvent(pActorSet, actor, time, LEAVING));
}


/**
 * Returns the set of exogenous events for the given period.
 */
const EventSet * Data::pEventSet(int period) const
{
	return this->lexogenousEvents[period];
}


// ----------------------------------------------------------------------------
// Section: Implementation of helper classes
// ----------------------------------------------------------------------------

/**
 * Returns if the first event is earlier than the second. If both events
 * happen simultaneously, the tie is broken by comparing the pointer values.
 */
bool EventComparator::operator()(const ExogenousEvent * pFirstEvent,
	const ExogenousEvent * pSecondEvent) const
{
	bool rc = false;

	if (pFirstEvent->time() != pSecondEvent->time())
	{
		rc = pFirstEvent->time() < pSecondEvent->time();
	}
	else
	{
		rc = pFirstEvent < pSecondEvent;
	}

	return rc;
}


// ----------------------------------------------------------------------------
// Section: Constraints
// ----------------------------------------------------------------------------

/**
 * Adds a new network constraint of the given type between two network
 * variables with the given names.
 */
const NetworkConstraint * Data::addNetworkConstraint(string networkName1,
	string networkName2,
	NetworkConstraintType type)
{
	NetworkConstraint * pConstraint =
		new NetworkConstraint(networkName1, networkName2, type);
	this->lnetworkConstraints.push_back(pConstraint);
	return pConstraint;
}


/**
 * Returns the vector of network constraints.
 */
const vector<const NetworkConstraint *> & Data::rNetworkConstraints() const
{
	return this->lnetworkConstraints;
}

}
