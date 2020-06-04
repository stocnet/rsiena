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
 * @param[in] parameter the value of the internal effect parameter this
 * function is associated with
 * @param[in] excludeMissing: whether to exclude missing values
 */
CovariateDistance2InAlterNetworkFunction::
CovariateDistance2InAlterNetworkFunction(string networkName, string
	covariateName, bool excludeMissing, bool total) :
	CovariateDistance2NetworkFunction(networkName, covariateName, excludeMissing, false)
{
	this->lexcludeMissing = excludeMissing;
	this->ltotal = total;
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double CovariateDistance2InAlterNetworkFunction::value(int alter)
{
	double value = 0;
	if (!(this->lexcludeMissing && this->missingInDummy(alter)))
	{
		if (ltotal)
		{
			value = this->totalInAlterValue(alter);
		}
		else
		{
			value = this->averageInAlterValue(alter);
		}
		int tieValue =  this->pNetwork()->tieValue(this->ego(), alter);
		if (tieValue == 1)
		{
			if (ltotal)
			{
				value = (value -
					CovariateNetworkAlterFunction::value(this->ego()));
			}
			else
			{
				int degree = this->pNetwork()->inDegree(alter);
				if (degree > 1)
				{
					value = (degree * value -
				CovariateNetworkAlterFunction::value(this->ego()))/(degree - 1);
				}
				else
				{
					value = CovariateNetworkAlterFunction::covmean();
				}
			}
		}
	}
	return value;
}


}
