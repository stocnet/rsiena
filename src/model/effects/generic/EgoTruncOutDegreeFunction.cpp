/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File:EgoTruncOutDegreeFunction.cpp
 *
 * Description: This file contains the implementation of the class
 *EgoTruncOutDegreeFunction.
 *****************************************************************************/

#include <cmath>
#include "EgoTruncOutDegreeFunction.h"
#include "data/NetworkLongitudinalData.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
EgoTruncOutDegreeFunction::EgoTruncOutDegreeFunction(string networkName, bool root, 
								bool threshold, int p) :
	NetworkAlterFunction(networkName)
{
	this->lroot = root;
	this->lthreshold = threshold;
	this->lsqrtTable = SqrtTable::instance();
	this->lintp = p;
	if (this->lroot)
	{
		this->lp = this->lsqrtTable->sqrt(p);
	}
	else
	{
		this->lp = p;
	}
	this->lvariableName = networkName;
// centering and root cannot occur simultaneously
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double EgoTruncOutDegreeFunction::value(int alter)  const
{
	double statistic = 0;
	double degree = this->pNetwork()->outDegree(this->ego());
	if (degree > this->lintp)
	{
		if (this->lthreshold)
		{
			statistic = 1;
		}
		else
		{
		if (this->lroot)
			{
				degree = this->lsqrtTable->sqrt(degree);
			}
			statistic = degree - this->lp;
		}
	}
	
	return statistic;
}

}
