#include <stdexcept>

#include "State.h"
#include "data/Data.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ContinuousLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/ContinuousVariable.h"
#include "model/settings/Setting.h"
#include "model/settings/PrimarySetting.h"
#include "network/OneModeNetwork.h"

using namespace std;

namespace siena
{

/**
 * Creates a state of variables as of the given observation of the given
 * Data object. The values may be copied or referenced directly as indicated
 * by the parameter <code>ownedValues</code>.
 */
State::State(const Data * pData, int observation, bool ownedValues)
{
	const vector<LongitudinalData *> & rVariables =
		pData->rDependentVariableData();

	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);
		ContinuousLongitudinalData * pContinuousData =
			dynamic_cast<ContinuousLongitudinalData *>(rVariables[i]);

		if (pNetworkData)
		{
			const Network * pNetwork = pNetworkData->pNetwork(observation);

			if (ownedValues)
			{
				pNetwork = pNetwork->clone();
			}

			this->lnetworks[pNetworkData->name()] = pNetwork;
		}
		else if (pBehaviorData)
		{
			const int * values = pBehaviorData->values(observation);

			if (ownedValues)
			{
				int * copies = new int[pBehaviorData->n()];

				for (int actor = 0; actor < pBehaviorData->n(); actor++)
				{
					copies[actor] = values[actor];
				}

				values = copies;
			}

			this->lbehaviors[pBehaviorData->name()] = values;
		}
		else if (pContinuousData)
		{
			const double * values = pContinuousData->values(observation);

			if (ownedValues)
			{
				double * copies = new double[pContinuousData->n()];

				for (int actor = 0; actor < pContinuousData->n(); actor++)
				{
					copies[actor] = values[actor];
				}

				values = copies;
			}

			this->lcontinuous[pContinuousData->name()] = values;
		}
		else
		{
			throw domain_error("unexpected class for longitudinal data: " + rVariables[i]->name());
		}
	}

	// SimVariables!
	const vector<LongitudinalData *> & rSimVariables = pData->rSimVariableData();
	for (unsigned i = 0; i < rSimVariables.size(); i++) {
		NetworkLongitudinalData * pNetworkData = dynamic_cast<NetworkLongitudinalData *>(rSimVariables[i]);
		if (pNetworkData) {
			const Network * pNetwork = pNetworkData->pNetwork(observation);
			if (ownedValues) {
				pNetwork = pNetwork->clone();
			}
			this->lnetworks[pNetworkData->name()] = pNetwork;
		} else {
			throw domain_error("unexpected class for simulated data: " + rSimVariables[i]->name());
		}
	}

	this->lownedValues = ownedValues;
}


State::State(EpochSimulation * pSimulation)
{
	const vector<DependentVariable *> & rVariables = pSimulation->rVariables();
	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkVariable * pNetworkVariable =
			dynamic_cast<NetworkVariable *>(rVariables[i]);
		BehaviorVariable * pBehaviorVariable =
			dynamic_cast<BehaviorVariable *>(rVariables[i]);

		if (pNetworkVariable)
		{
			this->lnetworks[pNetworkVariable->name()] =
				pNetworkVariable->pNetwork();

			const Setting * pSetting = pNetworkVariable->setting(1); // 0=universal, 1=primary
			if (pSetting) {
				const PrimarySetting * pPSetting = dynamic_cast<const PrimarySetting *>(pSetting);
				if (pPSetting) {
					if (pPSetting->pPrimaryNetwork() == 0) throw domain_error("no setting");
					std::string pname = "primary(" + pNetworkVariable->name() + ")";
					this->lnetworks[pname] = pPSetting->pPrimaryNetwork();
				}
			}
		}
		else if (pBehaviorVariable)
		{
			this->lbehaviors[pBehaviorVariable->name()] =
				pBehaviorVariable->values();
		}
		else
		{
			throw domain_error("unexpected class for dependent variable: " + rVariables[i]->name());
		}
	}

	const vector<ContinuousVariable *> & rContinuousVariables =
		pSimulation->rContinuousVariables();

	for (unsigned i = 0; i < rContinuousVariables.size(); i++)
	{
		this->lcontinuous[rContinuousVariables[i]->name()] =
				rContinuousVariables[i]->values();
	}

	this->lownedValues = false;
}


/**
 * Default constructor creating an empty state. The current values of dependent
 * variables can be stored later with the appropriate setters.
 */
State::State()
{
	this->lownedValues = false; // depends on the passed pointers
}


/**
 * Deallocates this state.
 */
State::~State()
{
	if (this->lownedValues)
	{
		this->deleteValues();
	}
}


/**
 * Returns the network for the given name or 0 if there is no such a network
 * stored in this state.
 */
const Network * State::pNetwork(string name) const
{
	map<string, const Network *>::const_iterator iter =
		this->lnetworks.find(name);
	if (iter != this->lnetworks.end())
	{
		return iter->second;
	}
	return 0;
}


/**
 * Stores the network for the given name.
 */
void State::pNetwork(string name, const Network * pNetwork)
{
	this->lnetworks[name] = pNetwork;
}


/**
 * Returns the values of the behavior variable with the given name, or 0
 * if no such values are stored in this state.
 */
const int * State::behaviorValues(string name) const
{
	const int * values = 0;
	map<string, const int *>::const_iterator iter =
		this->lbehaviors.find(name);

	if (iter != this->lbehaviors.end())
	{
		values = iter->second;
	}

	return values;
}


/**
 * Stores the values of a behavior variable with the given name.
 */
void State::behaviorValues(string name, const int * values)
{
	this->lbehaviors[name] = values;
}


/**
 * Returns the values of the continuous behavior variable with the given
 * name, or 0 if no such values are stored in this state.
 */
const double * State::continuousValues(string name) const
{
	const double * values = 0;
	map<string, const double *>::const_iterator iter =
		this->lcontinuous.find(name);

	if (iter != this->lcontinuous.end())
	{
		values = iter->second;
	}

	return values;
}


/**
 * Stores the values of a continuous behavior variable with the given name.
 */
void State::continuousValues(string name, const double * values)
{
	this->lcontinuous[name] = values;
}


/**
 * Deletes the values stored in this state (only called if lownedValues).
 */
void State::deleteValues()
{
	// Cannot use clearMap as the keys are not pointers.

	while (!this->lnetworks.empty())
	{
		const Network * pNetwork = this->lnetworks.begin()->second;
		this->lnetworks.erase(this->lnetworks.begin());
		delete pNetwork;
	}

	// Cannot use clearMap as the values are arrays and should be deleted with
	// delete[] operator.

	while (!this->lbehaviors.empty())
	{
		const int * values = this->lbehaviors.begin()->second;
		this->lbehaviors.erase(this->lbehaviors.begin());
		delete[] values;
	}

	while (!this->lcontinuous.empty())
	{
		const double * values = this->lcontinuous.begin()->second;
		this->lcontinuous.erase(this->lcontinuous.begin());
		delete[] values;
	}
}

}
