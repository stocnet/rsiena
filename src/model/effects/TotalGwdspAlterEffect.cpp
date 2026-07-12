/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalGwdspAlterEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwdspEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "TotalGwdspAlterEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/State.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include "model/tables/ConfigurationTable.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "NetworkDependentBehaviorEffect.h"
#include "model/EffectInfo.h"

// would probably be nice if GwespFunction could be used -> implement as network method?

using namespace std;

namespace siena
{


/**
 * Constructor.
 */
TotalGwdspAlterEffect::TotalGwdspAlterEffect(const EffectInfo * pEffectInfo, bool forward, bool nc) :
	NetworkDependentBehaviorEffect(pEffectInfo)
{
	this->lparameter = pEffectInfo->internalEffectParameter();
	this->lweight = -0.01 * this->lparameter;
	this->lexpmweight = exp(-this->lweight);
	this->lexpfactor = (1 - exp(this->lweight));
	this->lforward = forward;
	if (this->lparameter < 0)
	{
		throw runtime_error("Gwdsp must have nonnegative internal effect parameter");
	}
	this->lnc = nc;
}



/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void TotalGwdspAlterEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkDependentBehaviorEffect::initialize(pData, pState, period, pCache);
	// initialize the vector with weights for the GWDSP statistic
	// this is done several times during one estimation run (not elegant,
	// but not computationally burdensome either)
	double pow = 1;
	int m = this->pNetwork()->m();
	this->lcumulativeWeight.resize(m); // default values 0
	for (int i = 1; i < m; i++)
	{
		pow *= this->lexpfactor;
		this->lcumulativeWeight[i] = this->lexpmweight * (1 - pow);
	}

	this->lTwoPathCount.assign(this->pNetwork()->n(), 0);
	this->lTouched.clear();
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TotalGwdspAlterEffect::calculateChangeContribution(int actor,
		int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();
	
	if (pNetwork->outDegree(actor) > 0) 
	{
		this->lpInitialisedTable = 0;

		if (lforward)
					this->lpInitialisedTable = this->pTwoPathTable();
				else
					this->lpInitialisedTable = this->pInStarTable();

		double sumAlterValue = 0;
		for (int j = 0; j < this->n(); j++) //inefficient?
		{
			double alterValue = 0;
			int twoc = 0;
			if (j != actor)
			{
				if (this->lforward)
					twoc = this->lpInitialisedTable->get(j);
				else 
					twoc = this->lpInitialisedTable->get(j);
				alterValue = (this->lnc ? this->value(j) : 
										  this->centeredValue(j)) *
						this->lcumulativeWeight[twoc];
				sumAlterValue += alterValue;
			}
		}
		contribution = difference * sumAlterValue;
	}
	return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego.
 */
double TotalGwdspAlterEffect::egoStatistic(int ego, double * currentValues)
{
	const Network * pNetwork = this->pNetwork();
	// Scatter over ego's gateways h (out-ties) to accumulate pathCount(ego, j) for
	// each focal distance-2 actor j (FF: out-ties of h; FB: in-ties of h).
	for (IncidentTieIterator iterH = pNetwork->outTies(ego); iterH.valid(); iterH.next())
	{
		int h = iterH.actor();
		IncidentTieIterator iterJ = this->lforward ?
			pNetwork->outTies(h) : pNetwork->inTies(h);
		for (; iterJ.valid(); iterJ.next())
		{
			int j = iterJ.actor();
			if (j == ego) continue;
			if (this->lTwoPathCount[j] == 0) this->lTouched.push_back(j);
			this->lTwoPathCount[j]++;
		}
	}
	double statistic = 0;
	for (int j : this->lTouched)
	{
		statistic += (this->lnc ? currentValues[j] + this->overallCenterMean() :
								  currentValues[j]) *
				this->lcumulativeWeight[this->lTwoPathCount[j]];
		this->lTwoPathCount[j] = 0;
	}
	this->lTouched.clear();
	statistic *= this->lnc ?  currentValues[ego] + this->overallCenterMean() :
							  currentValues[ego];

	return statistic;
}

}
