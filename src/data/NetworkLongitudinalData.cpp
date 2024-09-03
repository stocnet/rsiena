/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkLongitudinalData.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkLongitudinalData class.
 *****************************************************************************/

#include <limits>
#include "NetworkLongitudinalData.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/NetworkUtils.h"
#include "network/IncidentTieIterator.h"
#include "data/ActorSet.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors, destructors.
// ----------------------------------------------------------------------------

/**
 * Creates a data object for storing the given number of observations of a
 * network. Initially the networks are empty at each observation.
 * @param[in] id the ID that is unique among all longitudinal data object
 * of the parent Data instance
 * @param[in] name the name of the corresponding network variable
 * @param[in] pSenders the set of actors acting as senders of ties
 * @param[in] pReceivers the set of actors acting as receivers of ties
 * @param[in] observationCount the number of observations to be stored
 */
NetworkLongitudinalData::NetworkLongitudinalData(int id,
	std::string name,
	const ActorSet * pSenders,
	const ActorSet * pReceivers,
	int observationCount,
	bool oneMode) :
		LongitudinalData(id, name, pSenders, observationCount)
{
	this->lpReceivers = pReceivers;
	this->lnetworks = new Network * [observationCount];
	this->lstructuralTieNetworks = new Network * [observationCount];
	this->lmissingTieNetworks = new Network * [observationCount];
	this->lnetworksLessMissings = new Network * [observationCount];
	this->lnetworksLessMissingStarts = new Network * [observationCount];
	this->lmaxDegree = std::numeric_limits<int>::max();
	this->lmodelType = 1;
	this->luniversalOffset = 0;
	this->ldensity = new double[observationCount];
	this->loneMode = oneMode;

	for (int i = 0; i < observationCount; i++)
	{
		if (oneMode)
		{
			this->lnetworks[i] = new OneModeNetwork(pSenders->n(), false);
			this->lstructuralTieNetworks[i] =
				new OneModeNetwork(pSenders->n(), false);
			this->lmissingTieNetworks[i] =
				new OneModeNetwork(pSenders->n(), false);
		}
		else
		{
			this->lnetworks[i] = new Network(pSenders->n(), pReceivers->n());
			this->lstructuralTieNetworks[i] =
				new Network(pSenders->n(), pReceivers->n());
			this->lmissingTieNetworks[i] =
				new Network(pSenders->n(), pReceivers->n());
		}
	}
}


/**
 * Deallocates this data object including all its networks.
 */
NetworkLongitudinalData::~NetworkLongitudinalData()
{
	for (int i = 0; i < this->observationCount(); i++)
	{
		delete this->lnetworks[i];
		delete this->lstructuralTieNetworks[i];
		delete this->lmissingTieNetworks[i];
		delete this->lnetworksLessMissings[i];
		delete this->lnetworksLessMissingStarts[i];
	}

	delete[] this->lnetworks;
	delete[] this->lstructuralTieNetworks;
	delete[] this->lmissingTieNetworks;
	delete[] this->ldensity;
	delete[] this->lnetworksLessMissings;
	delete[] this->lnetworksLessMissingStarts;

	this->lnetworks = 0;
	this->lstructuralTieNetworks = 0;
	this->lmissingTieNetworks = 0;
	this->ldensity = 0;
	this->lnetworksLessMissings = 0;
	this->lnetworksLessMissingStarts = 0;
}


// ----------------------------------------------------------------------------
// Section: Preprocessing
// ----------------------------------------------------------------------------

/**
 * Calculates various statistical properties from the stored network data. Also
 * store various versions of networks less missing or structurals for later use.
 */
void NetworkLongitudinalData::calculateProperties()
{
	// Calculate the overall average indegree and outdegree and
	// the observed network density at each observation.

	this->laverageInDegree = 0;
	this->laverageOutDegree = 0;
	this->laverageReciprocalDegree = 0;
	this->laverageSquaredInDegree = 0;
	this->laverageSquaredOutDegree = 0;

	for (int observation = 0;
		observation < this->observationCount();
		observation++)
	{
		Network * pNetwork = this->lnetworks[observation];
		Network * pMissingNetwork = this->lmissingTieNetworks[observation];

		for (int i = 0; i < this->lpReceivers->n(); i++)
		{
			this->laverageInDegree += pNetwork->inDegree(i);
			this->laverageSquaredInDegree +=
					(pNetwork->inDegree(i))*(pNetwork->inDegree(i));
		}

		int observedTieCount = 0;

		for (int i = 0; i < this->pActorSet()->n(); i++)
		{
			this->laverageOutDegree += pNetwork->outDegree(i);
			this->laverageSquaredOutDegree +=
				(pNetwork->outDegree(i))*(pNetwork->outDegree(i));
			observedTieCount +=
				pNetwork->outDegree(i) -
					commonActorCount(pNetwork->outTies(i),
						pMissingNetwork->outTies(i));
		}

		// Get the number of non-missing tie variables

		int nonMissingCount = this->n() * this->lpReceivers->n();

		if (this->loneMode)
		{
			const OneModeNetwork * pONetwork =
					dynamic_cast<const OneModeNetwork *>(pNetwork);
			if (!pONetwork)
			{
				throw logic_error("One-mode network expected in NetworkLongitudinalData.");
			}
			for (int i = 0; i < this->pActorSet()->n(); i++)
			{
				this->laverageReciprocalDegree += pONetwork->reciprocalDegree(i);
			}
		}
		else
		{
			// Don't count the diagonal entries for one-mode networks
			nonMissingCount -= this->n();
		}

		nonMissingCount -= pMissingNetwork->tieCount();

		// Calculate the observed density

		this->ldensity[observation] = 0;

		if (nonMissingCount > 0)
		{
			this->ldensity[observation] =
				((double) observedTieCount) / nonMissingCount;
		}
	}

	this->laverageInDegree /=
		this->lpReceivers->n() * this->observationCount();
	this->laverageOutDegree /=
		this->pActorSet()->n() * this->observationCount();
	this->laverageSquaredInDegree /=
		this->lpReceivers->n() * this->observationCount();
	this->laverageSquaredOutDegree /=
		this->pActorSet()->n() * this->observationCount();
	this->laverageReciprocalDegree /=
		this->pActorSet()->n() * this->observationCount();

	// data-less-missing-values is used in calculating statistics. Since it
	// does not change we store it, hoping we have enough space!
	for (int i = 0; i < this->observationCount(); i++)
	{
		this->lnetworksLessMissings[i] = this->pNetwork(i)->clone();
		this->lnetworksLessMissingStarts[i] = this->pNetwork(i)->clone();
		subtractNetwork(this->lnetworksLessMissings[i],
			this->pMissingTieNetwork(i));
		subtractNetwork(this->lnetworksLessMissingStarts[i],
			this->pMissingTieNetwork(i));
	}
	for (int i = 1; i < this->observationCount(); i++)
	{
		subtractNetwork(this->lnetworksLessMissings[i - 1],
			this->pMissingTieNetwork(i));
	}
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the set of actors acting as tie senders.
 */
const ActorSet * NetworkLongitudinalData::pSenders() const
{
	return this->pActorSet();
}


/**
 * Returns the set of actors acting as tie receivers.
 */
const ActorSet * NetworkLongitudinalData::pReceivers() const
{
	return this->lpReceivers;
}


/**
 * Returns the observed network as of the given observation moment.
 */
const Network * NetworkLongitudinalData::pNetwork(int observation) const
{
	return this->lnetworks[observation];
}


/**
 * Returns the network storing the structural tie indicators for the given
 * observation.
 */
const Network * NetworkLongitudinalData::pStructuralTieNetwork(int observation)
	const
{
	return this->lstructuralTieNetworks[observation];
}


/**
 * Returns the network storing the missing tie indicators for the given
 * observation.
 */
const Network * NetworkLongitudinalData::pMissingTieNetwork(int observation)
	const
{
	return this->lmissingTieNetworks[observation];
}

/**
 * Returns the network with missing values start or end zeroed, for the given
 * observation.
 */
const Network * NetworkLongitudinalData::pNetworkLessMissing(int observation)
	const
{
	return this->lnetworksLessMissings[observation];
}


/**
 * Returns the network with missing values at the start zeroed, for the given
 * observation.
 */
const Network * NetworkLongitudinalData::pNetworkLessMissingStart(int observation)
	const
{
	return this->lnetworksLessMissingStarts[observation];
}


/**
 * Returns the observed value of the tie from <i>i</i> to <i>j</i> at the given
 * observation.
 */
int NetworkLongitudinalData::tieValue(int i, int j, int observation) const
{
	return this->lnetworks[observation]->tieValue(i, j);
}


/**
 * Stores the observed value of the tie from <i>i</i> to <i>j</i> at the given
 * observation.
 */
void NetworkLongitudinalData::tieValue(int i,
	int j,
	int observation,
	int value)
{
	this->lnetworks[observation]->setTieValue(i, j, value);
}


/**
 * Returns if the tie value between the given actors is missing at the
 * given observation.
 */
bool NetworkLongitudinalData::missing(int i, int j, int observation) const
{
	return this->lmissingTieNetworks[observation]->tieValue(i, j);
}


/**
 * Stores if the tie value between the given actors is missing at the
 * given observation.
 */
void NetworkLongitudinalData::missing(int i, int j, int observation, bool flag)
{
	if (flag)
	{
		this->lmissingTieNetworks[observation]->setTieValue(i, j, 1);
	}
	else
	{
		this->lmissingTieNetworks[observation]->setTieValue(i, j, 0);
	}
}


/**
 * Returns if the tie value between the given actors is structurally determined
 * at the given observation.
 */
bool NetworkLongitudinalData::structural(int i, int j, int observation) const
{
	return this->lstructuralTieNetworks[observation]->tieValue(i, j);
}


/**
 * Stores if the tie value between the given actors is structurally determined
 * at the given observation.
 */
void NetworkLongitudinalData::structural(int i,
	int j,
	int observation,
	bool flag)
{
	if (flag)
	{
		this->lstructuralTieNetworks[observation]->setTieValue(i, j, 1);
	}
	else
	{
		this->lstructuralTieNetworks[observation]->setTieValue(i, j, 0);
	}
}


/**
 * Returns the number of structurally determined tie variables from the given
 * actor at the given observation.
 */
int NetworkLongitudinalData::structuralTieCount(int actor, int observation)
	const
{
	return this->lstructuralTieNetworks[observation]->outDegree(actor);
}

/**
 * Returns the number of receivers.
 */
int NetworkLongitudinalData::m() const
{
		if (this->loneMode )
		{
			return pActorSet()->n();
		}
		else
		{
			return this->lpReceivers->n();
		}
}

/**
 * Stores the maximum permitted out-degree of an actor.
 */
void NetworkLongitudinalData::maxDegree(int degree)
{
	this->lmaxDegree = degree;
}

/**
 * Stores the offset for the universal setting.
 */
void NetworkLongitudinalData::universalOffset(double offset) {
	this->luniversalOffset = offset;
}


/**
 * Stores the model type.
 */
void NetworkLongitudinalData::modelType(int type)
{
	this->lmodelType = type;
}

/**
 * Returns the maximum permitted out-degree of an actor.
 */
int NetworkLongitudinalData::maxDegree() const
{
	return this->lmaxDegree;
}

/**
 * Returns the offset for the universal setting.
 */
double NetworkLongitudinalData::universalOffset() const {
	return this->luniversalOffset;
}

/**
 * Returns the model type.
 */
int NetworkLongitudinalData::modelType() const
{
	return this->lmodelType;
}

/**
 * Returns whether the model type is NETCONTEMP.
 */
bool NetworkLongitudinalData::networkModelTypeContemp() const
{
	// This uses the order specified in 
	// enum NetworkModelType in DependentVariable.h
	return (this->lmodelType==11);
}
/**
 * Stores the average in-degree over all receivers and observations.
 */
void NetworkLongitudinalData::averageInDegree(double val)
{
	this->laverageInDegree = val;
}


/**
 * Store the average out-degree over all senders and observations.
 */
void NetworkLongitudinalData::averageOutDegree(double val)
{
	this->laverageOutDegree = val;
}

/**
 * Store the average reciprocal degree over all senders and observations.
 */
void NetworkLongitudinalData::averageReciprocalDegree(double val)
{
	this->laverageReciprocalDegree = val;
}

/**
 * Returns the average in-degree over all receivers and observations.
 */
double NetworkLongitudinalData::averageInDegree() const
{
	return this->laverageInDegree;
}

/**
 * Returns the average squared in-degree over all receivers and observations.
 */
double NetworkLongitudinalData::averageSquaredInDegree() const
{
	return this->laverageSquaredInDegree;
}


/**
 * Returns the average out-degree over all senders and observations.
 */
double NetworkLongitudinalData::averageOutDegree() const
{
	return this->laverageOutDegree;
}

/**
 * Returns the average squared out-degree over all senders and observations.
 */
double NetworkLongitudinalData::averageSquaredOutDegree() const
{
	return this->laverageSquaredOutDegree;
}

/**
 * Returns the average reciprocal degree over all senders and observations.
 * and 0 if two-mode
 */
double NetworkLongitudinalData::averageReciprocalDegree() const
{
	return this->laverageReciprocalDegree;
}

/**
 * Returns the relative frequency of the given value among the
 * observed value at the given observation.
 */
double NetworkLongitudinalData::observedDistribution(int value,
	int observation) const
{
	double frequency = 0;

	if (value == 1)
	{
		frequency = this->ldensity[observation];
	}
	else if (value == 0)
	{
		frequency = 1 - this->ldensity[observation];
	}

	return frequency;
}

/**
 * Returns whether this is a one mode network or not.
 */
bool NetworkLongitudinalData::oneModeNetwork() const
{
	return this->loneMode;
}
/**
 *Stores a setting name for this network.
 */
void NetworkLongitudinalData::addSettingName(const string id,
		const string settingType, const string covarName,
		const Permission_Type permType)
{
	this->lsettingNames.push_back(
			SettingInfo(id, settingType, covarName, permType));
}

/**
 * Returns the collection of settings names for this network.
 */
const std::vector<SettingInfo> & NetworkLongitudinalData::rSettingNames() const {
	return this->lsettingNames;
}
}
