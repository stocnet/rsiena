/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateAndNetworkBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateAndNetworkBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>
//#include "R_ext/Print.h"

#include "CovariateAndNetworkBehaviorEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/State.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
CovariateAndNetworkBehaviorEffect::CovariateAndNetworkBehaviorEffect(
	const EffectInfo * pEffectInfo) :
	CovariateDependentBehaviorEffect(pEffectInfo)
{
	// set up the extras if any
	this->laverageAlterValues = 0;
	this->lminimumAlterValues = 0;
	this->lmaximumAlterValues = 0;
	this->ltotalAlterValues = 0;
	this->laverageInAlterValues = 0;
	this->ltotalInAlterValues = 0;
	this->laverageAlterMissing = 0;
	this->laverageInAlterMissing = 0;
}

/**
 * Deallocates this effect object;
 */
CovariateAndNetworkBehaviorEffect::~CovariateAndNetworkBehaviorEffect()
{
	delete [] this->laverageAlterValues;
	delete [] this->lminimumAlterValues;
	delete [] this->lmaximumAlterValues;
	delete [] this->ltotalAlterValues;
	delete [] this->laverageInAlterValues;
	delete [] this->ltotalInAlterValues;
	delete [] this->laverageAlterMissing;
	delete [] this->laverageInAlterMissing;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateAndNetworkBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentBehaviorEffect::initialize(pData, pState, period, pCache);

	string networkName = this->pEffectInfo()->interactionName2();

	this->lpNetwork = pState->pNetwork(networkName);

	if (!this->lpNetwork)
	{
		throw logic_error("Network '" + networkName + "' expected.");
	}
	// initialize extras if any
	if (this->laverageAlterValues)
	{
		delete [] this->laverageAlterValues;
	}
	if (this->lminimumAlterValues)
	{
		delete [] this->lminimumAlterValues;
	}
	if (this->lmaximumAlterValues)
	{
		delete [] this->lmaximumAlterValues;
	}
	if (this->ltotalAlterValues)
	{
		delete [] this->ltotalAlterValues;
	}
	if (this->laverageInAlterValues)
	{
		delete [] this->laverageInAlterValues;
	}
	if (this->ltotalInAlterValues)
	{
		delete [] this->ltotalInAlterValues;
	}
	if (this->laverageAlterMissing)
	{
		delete[] this->laverageAlterMissing;
	}
	if (this->laverageInAlterMissing)
	{
		delete[] this->laverageInAlterMissing;
	}
	this->laverageAlterValues = new double[this->lpNetwork->n()];
	this->lminimumAlterValues = new double[this->lpNetwork->n()];
	this->lmaximumAlterValues = new double[this->lpNetwork->n()];
	this->ltotalAlterValues = new double[this->lpNetwork->n()];
	this->laverageInAlterValues = new double[this->lpNetwork->m()];
	this->ltotalInAlterValues = new double[this->lpNetwork->m()];
	this->laverageAlterMissing = new bool[this->lpNetwork->n()];
	this->laverageInAlterMissing = new bool[this->lpNetwork->m()];
}

/**
 * Returns if the dummy covariate value for the given actor is based on
 * all missing values.
 */
bool CovariateAndNetworkBehaviorEffect::missingDummy(int i) const
{
	return this->laverageAlterMissing[i];
}

/**
 * Returns if the dummy covariate value for the given actor is based on
 * all missing values.
 */
bool CovariateAndNetworkBehaviorEffect::missingInDummy(int i) const
{
	return this->laverageInAlterMissing[i];
}

/**
 * Returns the average alter covariate value for the given actor.
 */
double CovariateAndNetworkBehaviorEffect::averageAlterValue(int i) const
{
	return this->laverageAlterValues[i];
}

/**
 * Returns the minimum alter covariate value for the given actor.
 */
	double CovariateAndNetworkBehaviorEffect::minimumAlterValue(int i) const
	{
		return this->lminimumAlterValues[i];
	}

/**
 * Returns the maximum alter covariate value for the given actor.
 */
	double CovariateAndNetworkBehaviorEffect::maximumAlterValue(int i) const
	{
		return this->lmaximumAlterValues[i];
	}

/**
 * Returns the total alter covariate value for the given actor.
 */
double CovariateAndNetworkBehaviorEffect::totalAlterValue(int i) const
{
	return this->ltotalAlterValues[i];
}
/**
 * Returns the average in-alter covariate value for the given actor.
 */
double CovariateAndNetworkBehaviorEffect::averageInAlterValue(int i) const
{
	return this->laverageInAlterValues[i];
}

/**
 * Returns the total in-alter covariate value for the given actor.
 */
double CovariateAndNetworkBehaviorEffect::totalInAlterValue(int i) const
{
	return this->ltotalInAlterValues[i];
}

/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void CovariateAndNetworkBehaviorEffect::preprocessEgo(int ego)
{
	//CovariateDependentBehaviorEffect::preprocessEgo(ego);

	// set up the covariate based on current values of the network
	const Network * pNetwork = this->pNetwork();

	for (int i = 0; i < pNetwork->n(); i++)
	{
		this->laverageAlterMissing[i] = false;
		this->ltotalAlterValues[i] = 0;

		int counter = 1;
		if (pNetwork->outDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();

				if (!this->missingCovariate(j, this->period()))
				{
					if (counter == 1)
					{
						this->lminimumAlterValues[i] = this->covariateValue(j);
						this->lmaximumAlterValues[i] = this->covariateValue(j);
						counter++;
					}
					else
					{
						if (this->lminimumAlterValues[i] > this->covariateValue(j))
						{
							this->lminimumAlterValues[i] = this->covariateValue(j);
						}
						if (this->lmaximumAlterValues[i] < this->covariateValue(j))
						{
							this->lmaximumAlterValues[i] = this->covariateValue(j);
				}
					}
				}
				this->ltotalAlterValues[i] += this->covariateValue(j);
// 				Rprintf("%d %f %d %d %d %d\n",
// 					j,
// 					this->covariateValue(j),
// 					this->period(),
// 					this->missingCovariate(j, this->period()),
// 					numberNonMissing, i);
			}

			if(counter == 1)
			{
				this->lminimumAlterValues[i] = this->covariateMean();
				this->lmaximumAlterValues[i] = this->covariateMean();
				this->laverageAlterMissing[i] = true;
			}

			this->laverageAlterValues[i] =
					(this->ltotalAlterValues[i] / pNetwork->outDegree(i));
		}
		else
		{
			this->laverageAlterValues[i] = this->covariateMean();
			this->lminimumAlterValues[i] = this->covariateMean();
			this->lmaximumAlterValues[i] = this->covariateMean();
			this->ltotalAlterValues[i] = 0;
		}
//		Rprintf("%d %f\n", i,this->laverageAlterValues[i]);
	}

	for (int i = 0; i < pNetwork->m(); i++)
	{
		this->laverageInAlterMissing[i] = false;
		int numberNonMissing = 0;
		this->ltotalInAlterValues[i] = 0;
		if (pNetwork->inDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->inTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->ltotalInAlterValues[i] += this->covariateValue(j);
				if (!this->missingCovariate(j, this->period()))
				{
					numberNonMissing++;
				}
			}
			this->laverageInAlterValues[i] =
					(this->ltotalInAlterValues[i] / pNetwork->inDegree(i));
			if (numberNonMissing == 0)
			{
				this->laverageInAlterMissing[i] = true;
			}
		}
		else
		{
			this->laverageInAlterValues[i] = this->covariateMean();
			this->ltotalInAlterValues[i] = 0;
		}
	}
}

}
