/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * ConstantFunction.
 *****************************************************************************/

#include <stdexcept>
#include "ConstantFunction.h"
#include "data/NetworkLongitudinalData.h"
#include "data/Data.h"

using namespace std;

namespace siena
{

ConstantFunction::ConstantFunction(double constant)
{
	this->lconstant = constant;
	this->lconstantType = VALUE;
	this->lpFunction = 0;
}


ConstantFunction::ConstantFunction(string variableName,
	ConstantType constantType)
{
	this->lvariableName = variableName;
	this->lconstantType = constantType;
	this->lpFunction = 0;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ConstantFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);

	if (this->networkConstant())
	{
		NetworkLongitudinalData * pNetworkData =
			pData->pNetworkData(this->lvariableName);

		if (!pNetworkData)
		{
			throw logic_error(
				"Network data for " + this->lvariableName + " expected.");
		}

		if (this->lconstantType == AVERAGE_IN_DEGREE)
		{
			this->lconstant = pNetworkData->averageInDegree();
		}
		else if (this->lconstantType == AVERAGE_OUT_DEGREE)
		{
			this->lconstant = pNetworkData->averageOutDegree();
		}

		if (this->lpFunction)
		{
			this->lconstant = this->lpFunction(this->lconstant);
		}
	}
}


/**
 * Returns if this constant is to be read from the observed data of a network
 * variable.
 */
bool ConstantFunction::networkConstant() const
{
	return this->lconstantType == AVERAGE_IN_DEGREE ||
		this->lconstantType == AVERAGE_OUT_DEGREE;
}


double ConstantFunction::value(int alter)
{
	return this->lconstant;
}


void ConstantFunction::pFunction(double (* pFunction)(double))
{
	this->lpFunction = pFunction;
}

}
