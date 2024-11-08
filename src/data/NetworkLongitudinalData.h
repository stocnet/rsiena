/******************************************************************************
* SIENA: Simulation Investigation for Empirical Network Analysis
*
* Web: http://www.stats.ox.ac.uk/~snijders/siena/
*
* File: NetworkLongitudinalData.h
*
* Description: This file contains the definition of the
* NetworkLongitudinalData class.
*****************************************************************************/

#ifndef NETWORKLONGITUDINALDATA_H_
#define NETWORKLONGITUDINALDATA_H_

#include <vector>
#include "data/LongitudinalData.h"
#include "../model/settings/SettingInfo.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
* This class stores one or more observations of a network variable.
*/
class NetworkLongitudinalData : public LongitudinalData
{
public:
NetworkLongitudinalData(int id,
std::string name,
const ActorSet * pSenders,
const ActorSet * pReceivers,
int observationCount,
bool oneMode);
virtual ~NetworkLongitudinalData();

const ActorSet * pSenders() const;
const ActorSet * pReceivers() const;
const Network * pNetwork(int observation) const;
const Network * pStructuralTieNetwork(int observation) const;
const Network * pMissingTieNetwork(int observation) const;
const Network * pNetworkLessMissing(int observation) const;
const Network * pNetworkLessMissingStart(int observation) const;

int tieValue(int i, int j, int observation) const;
void tieValue(int i, int j, int observation, int value);
bool missing(int i, int j, int observation) const;
void missing(int i, int j, int observation, bool flag);
bool structural(int i, int j, int observation) const;
void structural(int i, int j, int observation, bool flag);
int structuralTieCount(int actor, int observation) const;

void calculateProperties();

int m() const; // number of receivers.
void maxDegree(int degree);
int maxDegree() const;
void universalOffset(double offset);
double universalOffset() const;
double averageInDegree() const; // average over all observations,
//calculated in calculateProperties
double averageOutDegree() const;// average over all observations,
//calculated in calculateProperties
double averageSquaredInDegree() const; // average over all observations,
//calculated in calculateProperties
double averageSquaredOutDegree() const;// average over all observations,
//calculated in calculateProperties
double averageReciprocalDegree() const;// average over all observations,
//calculated in calculateProperties; 0 if two-mode
void averageInDegree(double val);
void averageOutDegree(double val);
void averageReciprocalDegree(double val);
virtual double observedDistribution(int value, int observation) const;

void modelType(int type);
int modelType() const;
bool networkModelTypeContemp() const;

bool oneModeNetwork() const;

const std::vector<SettingInfo> & rSettingNames() const;
void addSettingName(const std::string id, const std::string settingType,
const std::string covarName, const Permission_Type permType);

private:
// The set of actors receiving the ties of the network
const ActorSet * lpReceivers;

// A collection of networks, one per each observation
Network ** lnetworks;

// Structural tie indicators per each observation stored as
// binary networks for efficiency

Network ** lstructuralTieNetworks;

// Missing tie indicators per each observation stored as
// binary networks for efficiency

Network ** lmissingTieNetworks;

// A collection of networks, one per each observation
Network ** lnetworksLessMissings;

// A collection of networks, one per each observation
Network ** lnetworksLessMissingStarts;

// The maximum permitted out-degree of an actor. Infinity by default.
int lmaxDegree {};

// The model type.
int lmodelType {};

// The offset for the non-primary settings. 0 by default.
double luniversalOffset {};

// The average in-degree over all receivers and observations
double laverageInDegree {};
double laverageSquaredInDegree {};

// The average out-degree (and squared out-degree) over all senders and observations
double laverageOutDegree {};
double laverageSquaredOutDegree {};

// The average reciprocal degree over all receivers and observations
double laverageReciprocalDegree {};

// The density of the network at each observation with the missing
// values excluded from the calculation.
double * ldensity {};

// Flag: is this data object bipartite or oneMode.
bool loneMode {};

// vector of setting names for this network
std::vector<SettingInfo> lsettingNames;
};

}
#endif /*NETWORKLONGITUDINALDATA_H_*/
