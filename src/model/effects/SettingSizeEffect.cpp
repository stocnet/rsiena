/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SettingSizeEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * SettingSizeEffect.
 *****************************************************************************/

#include <cmath>
#include "SettingSizeEffect.h"
#include "model/State.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{


/**
 * Constructor.
 */
SettingSizeEffect::SettingSizeEffect(const EffectInfo * pEffectInfo, bool difference,
		bool logar, bool root, bool inv, bool creation, bool evalDifference): SettingsNetworkEffect(pEffectInfo)
{
	this->lparameter = pEffectInfo->internalEffectParameter();
	this->ldifference = difference;
	this->llogar = logar;
	this->lroot = root;
	this->linv = inv;
	this->lcreation = creation;
	this->levalDifference = evalDifference;
	this->levalLog = (fabs(this->lparameter) < 1e-4);
	this->levalSqrt = (fabs(this->lparameter -2) < 1e-4);
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double SettingSizeEffect::calculateContribution(int alter) const
{
	double contribution = 0;
	if (this->lcreation)
	{
		if ((!this->outTieExists(alter)) && (this->stepType() == 1)) // primary setting
		{
			contribution = this->settingDegree() - this->outDegree();
			if (contribution >= 1)
			{
				if (this->llogar)
				{
					contribution = log(contribution+1);
				}
				else if (this->lroot)
				{
					contribution = sqrt(contribution);
				}
				else if (this->linv)
				{
					contribution = 1/(contribution + 1);
				}
			}
		}
	}
	else
	{
		contribution = this->settingDegree();
		if (this->ldifference)
		{
			contribution -= this->outDegree();
		}			
		if (contribution >= 1)
		{
			if (this->llogar)
			{
				contribution = log(contribution+1);
			}
			else if (this->lroot)
			{
				contribution = sqrt(contribution);
			}
			else if (this->linv)
			{
				contribution = 1/(contribution + 1);
			}
		}
	}
	return contribution;
}


/**
 * The contribution of ego to the statistic.
 */
double SettingSizeEffect::egoStatistic(int ego, const Network * pNetwork)
{
	double statistic = 0;
	int size = this->settingDegree();
	if (this->levalDifference)
	{
		size -= this->outDegree();
	}
	if (size > 0)
	{
		if (this->levalLog)
		{
			statistic = log(size);
		}
		else if (this->levalSqrt)
		{
			statistic = sqrt(size);
		}
		else if (this->linv)
		{
			statistic = 1/(size+1);
		}
		else
		{
			statistic = size;
		}
	}
	else
	{
		statistic = 0;
	}
	return statistic;
}

}
