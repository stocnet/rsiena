/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateInTiesFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * SameCovariateInTiesFunction.
 *****************************************************************************/

#include <cmath>
#include "SameCovariateInTiesFunction.h"
#include "NetworkAlterFunction.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"
#include "utils/Utils.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 * @param[in] covariateName the name of the covariate this function is
 * associated with
 * @param[in] excludeMissing whether to exclude missing values
 */
SameCovariateInTiesFunction::SameCovariateInTiesFunction(
		string networkName, string covariateName, bool sameValue,
		bool sameVariable, int parameter, bool excludeMissing) :
	CovariateNetworkAlterFunction(networkName, covariateName)
{
	this->lsameValue = sameValue;
	this->lsameVariable = sameVariable;
	this->lexcludeMissing = excludeMissing;
	this->laverage = (parameter >= 3);
	this->lroot = ((parameter == 2)||(parameter == 4));
	this->lCovNumberEgo = 1;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void SameCovariateInTiesFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateNetworkAlterFunction::initialize(pData, pState, period, pCache);	
	
	if (this->laverage)
	{
		int covMin = int(round(this->covariateMinimum()));
		int covMax = int(round(this->covariateMaximum()));
		covMax++;
	
		if (covMin < 0)
		{
			throw logic_error("sameXInPop: minimum of covariate is negative");		
		}
		if (covMax > 20)
		{
			throw logic_error("sameXInPop: covariate has a maximum which is too large");		
		}
		this->lpCovariateNumbers = new int[covMax] {};
		
		for (int i = 0; i < covMax; i++)
		{
			this->lpCovariateNumbers[i] = 0;		
		}
		for (int i = 0; i < this->covariateN(); i++)
		{
			this->lpCovariateNumbers[this->covIntValue(i)]++;		
		}
	}
}



/**
 * Deallocates this SameCovariateInTiesFunction object.
 */
SameCovariateInTiesFunction::~SameCovariateInTiesFunction()
{
	if (this->laverage)
	{
		delete[] this->lpCovariateNumbers;
		this->lpCovariateNumbers = 0;
	}
}

/**
 * Does the necessary preprocessing work for calculating the
 * predicate for a specific ego. This method must be invoked before
 * calling SameCovariateInTiesFunction::value(...).
 */
void SameCovariateInTiesFunction::preprocessEgo(int ego)
{
	CovariateNetworkAlterFunction::preprocessEgo(ego);
	this->lCovEgo = this->covvalue(ego);
// Not rounded to integer, because for !this->laverage,
// the function works also for non-integer covariates.
	if (this->laverage)
	{
		this->lCovNumberEgo = this->lpCovariateNumbers[this->covIntValue(ego)];
	}
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double SameCovariateInTiesFunction::value(int alter) const
{
	double statistic = 0;
	double vval = 0;
	if  (!(this->lexcludeMissing && this->missing(this->ego())))
	{
		const Network * pNetwork = this->pNetwork();
		// Iterate over incoming ties of alter
		if (this->lsameValue)
		{
			for (IncidentTieIterator iter =	pNetwork->inTies(alter);
			iter.valid();
			iter.next())
			{
				// Get the sender of the incoming tie.
				int h = iter.actor();
				// ego needs to have the same covariate value as h:
				if (!(this->lexcludeMissing && this->missing(h)))
				{
					if ((fabs(this->covvalue(h)
									- this->lCovEgo) < EPSILON))
					{
						statistic++;
					}
				}
			}
			// Add the following just like for inPop effect:
			if ((this->lsameVariable) && (!this->outTieExists(alter)))
			{
				// The statistic will increase after introducing the tie
				statistic++;
			}
		}
		else
		{
			for (IncidentTieIterator iter =	pNetwork->inTies(alter);
			iter.valid();
			iter.next())
			{
				// Get the sender of the incoming tie.
				int h = iter.actor();
				// ego needs to have a different covariate value than h:
				if (!(this->lexcludeMissing && this->missing(h)))
				{
					if (fabs(this->covvalue(h)
									- this->lCovEgo) >= EPSILON)
					{
						statistic++;
					}
				}
			}
		}
		if (this->laverage)
		{
			if (this->lsameValue)
			{
				vval = statistic/(this->lCovNumberEgo); 
			}
			else
			{
				int totdiff = this->covariateN() - this->lCovNumberEgo ;
				if (totdiff > 0) // else 0/0
				{
					vval = statistic/totdiff; 
				}
			}
		}
		else
		{
			vval = statistic;
		}
		if (this->lroot)
		{
			vval = sqrt(vval);
		}
	}
	return vval;
}

}
