/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2AlterNetworkFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateDistance2AlterNetworkFunction.
 *****************************************************************************/
#include <R_ext/Print.h>
#include "CovariateDistance2AlterNetworkFunction.h"
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
CovariateDistance2AlterNetworkFunction::
CovariateDistance2AlterNetworkFunction(string networkName, string
	covariateName, double parameter,  bool excludeMissing, bool total, bool nc) :
	CovariateDistance2NetworkFunction(networkName, covariateName, excludeMissing, true, nc)
{
	this->lparameter = parameter;
	this->lexcludeMissing = excludeMissing;
	this->ltotal = total;
}

/**
 * Deallocates this effect object;
 */
	CovariateDistance2AlterNetworkFunction::~CovariateDistance2AlterNetworkFunction()
{
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double CovariateDistance2AlterNetworkFunction::value(int alter) const
{
//	Rprintf("cccc %d %d\n", this->ego(), alter);
	double value = 0;
	if (!(this->lexcludeMissing && this->missingDummy(alter)))
	{
		if (ltotal)
		{
			value = this->totalAlterValue(alter);
		}
		else
		{
			value = this->averageAlterValue(alter);
		}
		if (this->lparameter == 2)
		{
			int tieValue =  this->pNetwork()->tieValue(alter, this->ego());
			if (tieValue == 1)
			{
				double egoVal = CovariateNetworkAlterFunction::resolvedCovvalue(this->ego(), this->nc());
				int degree = this->pNetwork()->outDegree(alter);
				//			Rprintf("before %d %f %d %f\n", degree, value,
				//this->ego(), CovariateDistance2NetworkFunction::value(this->ego()) );
				if (degree > 1)
				{
					if (ltotal)
					{
						value = (value - egoVal);
					}
					else
					{
						value = (degree * value - egoVal)/(degree - 1);
					}
				}
				else
				{
					value = this->nc() ? 0 : CovariateNetworkAlterFunction::covmean();
				}
				//Rprintf("stat after %d %f %d %f\n", degree, value,
				//	this->ego(),
				//	CovariateDistance2NetworkFunction::value(this->ego()) );
			}
		}
	}
//	Rprintf("stat r  %f f\n", value);

	return value;
}


}
