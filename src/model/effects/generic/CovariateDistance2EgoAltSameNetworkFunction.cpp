/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2EgoAltSameNetworkFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateDistance2EgoAltSameNetworkFunction.
 *****************************************************************************/
#include <R_ext/Print.h>
#include "CovariateDistance2EgoAltSameNetworkFunction.h"
#include "CovariateNetworkAlterFunction.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include <cmath>

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
CovariateDistance2EgoAltSameNetworkFunction::
CovariateDistance2EgoAltSameNetworkFunction(string networkName, string
	covariateName, bool excludeMissing, bool outgoing, double parameter) :
	CovariateDistance2NetworkFunction(networkName, covariateName, 
								excludeMissing, outgoing)
{
	this->lexcludeMissing = excludeMissing;
	this->loutgoing = outgoing;
	this->ltrunc = (std::round(parameter) == 0);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double CovariateDistance2EgoAltSameNetworkFunction::value(int alter) const
{
	double value = 0;
	IncidentTieIterator iter;
	int deg = 0;
	const Network * pNetwork = this->pNetwork();
	if (loutgoing)
	{
		iter = pNetwork->outTies(alter);
		deg = pNetwork->outDegree(alter);
	}
	else
	{
		iter = pNetwork->inTies(alter);		
		deg = pNetwork->inDegree(alter);
	}

	if ((deg > 0) && (!missing(this->ego())))
	{
		double egoValue = CovariateNetworkAlterFunction::covvalue(this->ego());
		int numberUsed = 0;
		for ( ; iter.valid(); iter.next())
		{
			int h = iter.actor();
			if ((!(this->lexcludeMissing && this->missing(h))) && 
												(!(h == this->ego())))
			{
				if (CovariateNetworkAlterFunction::covvalue(h) == egoValue)
				{
					value++;
				}
				numberUsed++;
			}
		}
		
		if (ltrunc)
		{
			if (value > 0)
			{
				value = 1;			
			}
		}
		else
		{			
			if (numberUsed > 0)
			{
				value/=numberUsed;			
			}
		}		
	}
	return value;
}
}
