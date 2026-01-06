/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAvSimEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AltersCovariateAvRecAltEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>

#include "AltersCovariateAvRecAltEffect.h"
#include "data/Data.h"
#include "network/OneModeNetwork.h"
#include "network/CommonNeighborIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AltersCovariateAvRecAltEffect::AltersCovariateAvRecAltEffect(
	const EffectInfo * pEffectInfo, bool divide, bool same) :
	CovariateAndNetworkBehaviorEffect(pEffectInfo)
{
	this->ldivide = divide;
	// ldivide indicates whether there will be division by the outdegree of ego
	this->lsame = same;
	// same indicates whether the influence will be from same-X alters
}

/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void AltersCovariateAvRecAltEffect::preprocessEgo(int ego)
{
// do not repeat all this:	CovariateAndNetworkBehaviorEffect::preprocessEgo(ego);
	this->lTotalAlterValue = 0;
	int neighborCount = 0;
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());
	if (!pNetwork)
	{
		throw runtime_error(string("One-mode network expected in ") +
			"***RecAltAlt***X effect");
	}

	if (pNetwork->reciprocalDegree(ego) > 0)
	{
		double covEgo = this->covariateValue(ego);

		for (CommonNeighborIterator iter = pNetwork->reciprocatedTies(ego);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();                // identifies alter
			if (this->lsame)
			{
				if (fabs(this->covariateValue(j) - covEgo) < 1e-6)
				{
					this->lTotalAlterValue += this->centeredValue(j);
					neighborCount++;	
				}
			}
			else
			{
				this->lTotalAlterValue += this->centeredValue(j);
				neighborCount++;	
			}			
		}

		if ((this->ldivide) & (neighborCount > 0))
		{
			this->lTotalAlterValue /= neighborCount;
		}
	}

}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AltersCovariateAvRecAltEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->lTotalAlterValue;
}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AltersCovariateAvRecAltEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());
	int neighborCount = 0;
	double covEgo = this->covariateValue(ego);
	bool misEgo = (this->missing(this->period() + 1, ego) || 				
						this->missing(this->period() + 1, ego));

	if (!(this->lsame && misEgo))
	{
		for (CommonNeighborIterator iter = pNetwork->reciprocatedTies(ego);
			iter.valid();
			iter.next())
		{
			int j = iter.actor();

			if (!this->missing(this->period(), j) &&
				!this->missing(this->period() + 1, j) &&
				!this->missingCovariate(j,this->period()))
			{
				if (this->lsame) 
				{
					if (fabs(this->covariateValue(j) - covEgo) < 1e-6)
					{
						statistic += currentValues[j];	
						neighborCount++;	
					}
				}
				else
				{
					statistic += currentValues[j];
					neighborCount++;	
				}			
			}
		}
		statistic *= currentValues[ego];
		if ((neighborCount > 0) && (this->ldivide))
		{
			statistic /= neighborCount;
		}
	}

	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AltersCovariateAvRecAltEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;
	const OneModeNetwork * pNetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());
	
	if ((difference[ego] > 0) && (pNetwork->reciprocalDegree(ego) > 0)) 
	{
		int neighborCount = 0;
		double covEgo = this->covariateValue(ego);
		bool misEgo = (this->missing(this->period() + 1, ego) || 				
						this->missing(this->period() + 1, ego));

		if (!(this->lsame && misEgo))
		{
			for (CommonNeighborIterator iter = pNetwork->reciprocatedTies(ego);
				iter.valid();
				iter.next())
			{
				int j = iter.actor();

				if (!this->missing(this->period(), j) &&
					!this->missing(this->period() + 1, j) &&
					!this->missingCovariate(j,this->period()))
				{
					if (this->lsame) 
					{
						if (fabs(this->covariateValue(j) - covEgo) < 1e-6)
						{
							statistic += currentValues[j];	
							neighborCount++;			
						}
					}
					else
					{
						statistic += currentValues[j];
						neighborCount++;
					}
				}
			}
			statistic *= difference[ego];
			if ((neighborCount > 0) && (this->ldivide))
			{
				statistic /= neighborCount;
			}
		}
	}
	return -statistic;
}


}
