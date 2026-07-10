/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2InAlterNetworkFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateDistance2InAlterNetworkFunction.
 *****************************************************************************/
#include <R_ext/Print.h>
#include "CovariateDistance2InAlterNetworkFunction.h"
#include "CovariateNetworkAlterFunction.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 * @param[in] covariateName the name of the covariate this function is
 * associated with
 * @param[in] excludeMissing: whether to exclude missing values
 * @param[in] total: whether to calculate the total value
 * @param[in] nc: whether to use raw values (non-centered) of a behavior
 * variable instead of centered values
 */
CovariateDistance2InAlterNetworkFunction::
CovariateDistance2InAlterNetworkFunction(string networkName, string
	covariateName, bool excludeMissing, bool total, bool nc) :
	CovariateDistance2NetworkFunction(networkName, covariateName, excludeMissing, false, nc)
{
	this->lexcludeMissing = excludeMissing;
	this->ltotal = total;
}



void CovariateDistance2InAlterNetworkFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDistance2NetworkFunction::initialize(pData, pState, period, pCache);
}

/**
 * Does the necessary preprocessing work for calculating the
 * predicate for a specific ego. This method must be invoked before
 * calling CovariateDistance2InAlterNetworkFunction::value(...).
 */
void CovariateDistance2InAlterNetworkFunction::preprocessEgo(int ego)
{
	CovariateDistance2NetworkFunction::preprocessEgo(ego);
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double CovariateDistance2InAlterNetworkFunction::value(int alter) const
{
	double value = 0;
	if (!(this->lexcludeMissing && this->missingInDummy(alter)))
	{
		if (this->ltotal)
		{
			value = this->totalInAlterValue(alter);
		}
		else
		{
			value = this->averageInAlterValue(alter);
		}
		
		if (this->pNetwork()->tieValue(this->ego(), alter) == 1)
		{
			double egoVal = this->nc()
				? CovariateNetworkAlterFunction::uncenteredCovvalue(this->ego())
				: CovariateNetworkAlterFunction::covvalue(this->ego());
			if (this->ltotal)
			{
				value = (value - egoVal);
			}
			else
			{
				int degree = this->pNetwork()->inDegree(alter);
				if (degree > 1)
				{
					value = (degree * value - egoVal) / (degree - 1);
				}
				else
				{
					value = this->nc() ? 0 : CovariateNetworkAlterFunction::covmean();
				}
			}
		}
	}
	return value;
}


}
