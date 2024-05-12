/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Data.h
 *
 * Description: This file contains the definition of the Data class.
 *****************************************************************************/

#ifndef DATA_H_
#define DATA_H_

#include <vector>
#include <set>
#include <map>
#include <string>
#include "NetworkConstraint.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class LongitudinalData;
class NetworkLongitudinalData;
class OneModeNetworkLongitudinalData;
class BehaviorLongitudinalData;
class ContinuousLongitudinalData;
class ConstantCovariate;
class ChangingCovariate;
class ConstantDyadicCovariate;
class ChangingDyadicCovariate;
class ExogenousEvent;
class ActorSet;


// ----------------------------------------------------------------------------
// Section: Helper definitions
// ----------------------------------------------------------------------------

/**
 * This class defines a comparator class for sorting exogenous events
 * chronologically.
 */
class EventComparator
{
public:
	bool operator()(const ExogenousEvent * pFirstEvent,
		const ExogenousEvent * pSecondEvent) const;
};


/**
 * Shorthand for sorted sets of exogenous events of composition change.
 */
typedef std::set<ExogenousEvent *, EventComparator> EventSet;


// ----------------------------------------------------------------------------
// Section: Data class
// ----------------------------------------------------------------------------

/**
 * This class acts as a storage of the whole data set subject to actor-based
 * modeling.
 */
class Data
{
public:
	Data(int observationCount);
	virtual ~Data();

	inline int observationCount() const;

	const ActorSet * createActorSet(std::string name, int n);
	NetworkLongitudinalData * createNetworkData(std::string name,
		const ActorSet * pSenders,
		const ActorSet * pReceivers);
	OneModeNetworkLongitudinalData * createOneModeNetworkData(
			std::string name, const ActorSet * pActors);
	OneModeNetworkLongitudinalData * createOneModeSimNetworkData(
			std::string name, const ActorSet * pActors);
	BehaviorLongitudinalData * createBehaviorData(std::string name,
		const ActorSet * pActorSet);
	ContinuousLongitudinalData * createContinuousData(std::string name,
		const ActorSet * pActorSet);
	ConstantCovariate * createConstantCovariate(std::string name,
		const ActorSet * pActorSet);
	ChangingCovariate * createChangingCovariate(std::string name,
		const ActorSet * pActorSet);
	ConstantDyadicCovariate * createConstantDyadicCovariate(std::string name,
		const ActorSet * pFirstActorSet,
		const ActorSet * pSecondActorSet);
	ChangingDyadicCovariate * createChangingDyadicCovariate(std::string name,
		const ActorSet * pFirstActorSet,
		const ActorSet * pSecondActorSet);

	const ActorSet * pActorSet(std::string name) const;
	NetworkLongitudinalData * pNetworkData(std::string name) const;
	NetworkLongitudinalData * pSimNetworkData(std::string name) const;
	OneModeNetworkLongitudinalData * pOneModeNetworkData(std::string name) const;
	OneModeNetworkLongitudinalData * pOneModeSimNetworkData(std::string name) const;
	BehaviorLongitudinalData * pBehaviorData(std::string name) const;
	ContinuousLongitudinalData * pContinuousData(std::string name) const;
	ConstantCovariate * pConstantCovariate(std::string name) const;
	ChangingCovariate * pChangingCovariate(std::string name) const;
	ConstantDyadicCovariate * pConstantDyadicCovariate(std::string name) const;
	ChangingDyadicCovariate * pChangingDyadicCovariate(std::string name) const;

	const std::vector<const ActorSet *> & rActorSets() const;
	const std::vector<LongitudinalData *> & rDependentVariableData() const;
	const std::vector<LongitudinalData *> & rSimVariableData() const;
	const std::vector<ConstantCovariate *> & rConstantCovariates() const;
	const std::vector<ChangingCovariate *> & rChangingCovariates() const;
	const std::vector<ConstantDyadicCovariate *> & rConstantDyadicCovariates() const;
	const std::vector<ChangingDyadicCovariate *> & rChangingDyadicCovariates() const;

	void active(const ActorSet * pActorSet,
		int actor,
		int observation,
		bool flag);
	bool active(const ActorSet * pActorSet, int actor, int observation);
	void addJoiningEvent(int period,
		const ActorSet * pActorSet,
		int actor,
		double time);
	void addLeavingEvent(int period,
		const ActorSet * pActorSet,
		int actor,
		double time);
	const EventSet * pEventSet(int period) const;

	// Network constraints

	const NetworkConstraint * addNetworkConstraint(std::string networkName1,
		std::string networkName2,
		NetworkConstraintType type);
	const std::vector<const NetworkConstraint *> & rNetworkConstraints() const;

private:
	// The number of observations
	int lobservationCount {};

	// A collection of actor sets
	std::vector<const ActorSet *> lactorSets;

	// A collection of longitudinal data objects for all kinds of
	// dependent variables

	std::vector<LongitudinalData *> ldependentVariableData;
	std::vector<LongitudinalData *> lsimVariableData;

	// A collection of longitudinal continuous dependent variables
	// Necessary for computation of statistics based on all cont DVs
	std::vector<ContinuousLongitudinalData *> lcontinuousDependentVariableData;

	// A collection of constant covariates
	std::vector<ConstantCovariate *> lconstantCovariates;

	// A collection of changing covariates
	std::vector<ChangingCovariate *> lchangingCovariates;

	// A collection of constant dyadic covariates
	std::vector<ConstantDyadicCovariate *> lconstantDyadicCovariates;

	// A collection of changing dyadic covariates
	std::vector<ChangingDyadicCovariate *> lchangingDyadicCovariates;

	// lactive[s][i][k] indicates if actor i of the actor set s is active at
	// the observation k

	std::map<const ActorSet *, bool **> lactive;

	// A sorted set of exogenous events per each period
	std::vector<EventSet *> lexogenousEvents;

	// Network constraints like higher(network1,network2),
	// disjoint(network1,network2), etc.

	std::vector<const NetworkConstraint *> lnetworkConstraints;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the number of observations stored in this data object.
 */
inline int Data::observationCount() const
{
	return this->lobservationCount;
}

}

#endif /*DATA_H_*/
