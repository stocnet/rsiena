/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DiffusionRateEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * DiffusionRateEffect.
 *****************************************************************************/

#include <cmath>
#include <cstring>
#include "DiffusionRateEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/DiffusionEffectValueTable.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/IncidentTieIterator.h"
#include <R_ext/Print.h>

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pVariable the network variable this effect depends on
 * @param[in] pBehaviorVariable the behavior variable this effect depends on
 * @param[in] parameter the statistical parameter of this effect
 */
DiffusionRateEffect::DiffusionRateEffect(const NetworkVariable * pVariable,
	const BehaviorVariable * pBehaviorVariable,
	string effectName,
	double parameter,
	double internalEffectParameter)
{
	this->lpVariable = pVariable;
	this->lpBehaviorVariable = pBehaviorVariable;
	this->leffectName = effectName;
	double possibleDegreeNumer = this->lpBehaviorVariable->range()
		* max(this->lpVariable->n(), this->lpVariable->m());
	double possibleDegreeDenom = 1;

	if (effectName == "avExposure")
	{
		possibleDegreeDenom = max(this->lpVariable->n(), this->lpVariable->m());
	}
	if (effectName == "susceptAvIn")
	{
		possibleDegreeNumer *= max(this->lpVariable->n(), this->lpVariable->m());
		possibleDegreeDenom = max(this->lpVariable->n(), this->lpVariable->m());
	}
	if ((effectName == "infectDeg") || (effectName == "infectIn") ||
										(effectName == "infectOut"))
	{
		possibleDegreeNumer *= max(this->lpVariable->n(),
				this->lpVariable->m());
	}
	this->lpTable = new DiffusionEffectValueTable(possibleDegreeNumer,
			possibleDegreeDenom);
	this->lpTable->parameter(parameter);
	this->linternalEffectParameter = round(internalEffectParameter);
	this->labsInternalEffectParameter = std::abs(this->linternalEffectParameter);
	this->linternalNonZero = (this->linternalEffectParameter != 0);

	if (((effectName == "infectDeg") || (effectName == "infectIn") ||
				(effectName == "infectOut")) && (this->linternalEffectParameter < 0))
	{
		throw logic_error("Negative internal parameter not permitted for effect "+effectName);
	}
}

DiffusionRateEffect::DiffusionRateEffect(const NetworkVariable * pVariable,
	const BehaviorVariable * pBehaviorVariable,
	const ConstantCovariate * pConstantCovariate,
	const ChangingCovariate * pChangingCovariate,
	string effectName,
	double parameter,
	double internalEffectParameter)
{
	this->lpVariable = pVariable;
	this->lpBehaviorVariable = pBehaviorVariable;
	this->lpChangingCovariate = pChangingCovariate;
	this->lpConstantCovariate = pConstantCovariate;
	this->leffectName = effectName;
	this->linternalEffectParameter = round(internalEffectParameter);
	this->labsInternalEffectParameter = std::abs(this->linternalEffectParameter);
	this->linternalNonZero = (this->linternalEffectParameter != 0);

	double possibleDegreeNumer = 1;
	double possibleDegreeDenom = 1;

	if (effectName == "susceptAvCovar")
	{
		possibleDegreeNumer = (this->lpBehaviorVariable->range())
			* max(this->lpVariable->n(), this->lpVariable->m());
		possibleDegreeDenom = max(this->lpVariable->n(), this->lpVariable->m());
	}

	this->lpTable = new DiffusionEffectValueTable(possibleDegreeNumer,
		possibleDegreeDenom);
	this->lpTable->parameter(parameter);

	if ((effectName == "infectCovar") && (this->linternalEffectParameter < 0))
	{
		throw logic_error("Negative internal parameter not permitted for effect "+effectName);
	}
}

/**
 * Destructor.
 */
DiffusionRateEffect::~DiffusionRateEffect()
{
	delete this->lpTable;
	this->lpTable = 0;
}

/**
 * Returns the value of a proximity measure effect, a type of diffusion rate
 * effect.
 */
double DiffusionRateEffect::proximityValue(Network * pNetwork, int i,
		int egoNumer, int egoDenom) const
{
	int totalAlterValue = 0;
	int numInfectedAlter = 0;
	if (pNetwork->outDegree(i) > 0)
	{
		for (IncidentTieIterator iter = pNetwork->outTies(i);
			 iter.valid();
			 iter.next())
		{
			double alterValue = this->lpBehaviorVariable->value(iter.actor());

			if (alterValue >= 0.5)
			{
				numInfectedAlter++;
			}

			if (this->leffectName == "infectIn")
			{
				alterValue *= pNetwork->inDegree(iter.actor());
			}
			else if ((this->leffectName == "infectDeg") ||
						(this->leffectName == "infectOut"))
			{
				alterValue *= pNetwork->outDegree(iter.actor());
			}

			totalAlterValue += alterValue;
		}
	}

	if (this->linternalNonZero)
	{
		if (numInfectedAlter < this->labsInternalEffectParameter)
		{
			totalAlterValue = 0;
		}
		else if (this->linternalEffectParameter < 0)
		{
			if (totalAlterValue > this->labsInternalEffectParameter)
			{
				totalAlterValue = this->labsInternalEffectParameter;
			}
		}
	}

	totalAlterValue *= egoNumer;

	if (totalAlterValue == 0)
	{
		return 1; // = exp(0)
	}
	else
	{
		return this->lpTable->value(totalAlterValue, egoDenom);
	}
}

/**
 * Returns the contribution of this effect for the given actor.
 */

double DiffusionRateEffect::value(int i, int period) const
{
	Network * pNetwork = this->lpVariable->pNetwork();

	if (this->leffectName == "avExposure")
	{
		return this->proximityValue(pNetwork, i, 1, max(1,
				pNetwork->outDegree(i)));
	}
	else if (this->leffectName == "susceptAvIn")
	{
		return this->proximityValue(pNetwork, i, pNetwork->inDegree(i),
			max(1, pNetwork->outDegree(i)));
	}
	else if (this->leffectName == "totExposure" ||
		this->leffectName == "infectDeg" ||
		this->leffectName == "infectIn" ||
		this->leffectName == "infectOut")
	{
		return this->proximityValue(pNetwork, i, 1, 1);
	}
	else if (this->leffectName == "susceptAvCovar")
	{
		if (this->lpConstantCovariate)
		{
			return pow(this->proximityValue(
						pNetwork, i, 1, max(1, pNetwork->outDegree(i))),
					this->lpConstantCovariate->value(i));
		}
		else if (this->lpChangingCovariate)
		{
			return pow(this->proximityValue(pNetwork, i, 1, max(1,
						pNetwork->outDegree(i))),
				this->lpChangingCovariate->value(i,period));
		}
		else
			throw logic_error(
				"No individual covariate found.");
	}
	else if (this->leffectName == "infectCovar")
	{
		double totalAlterValue = 0;
		if (pNetwork->outDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				double alterValue = this->lpBehaviorVariable->
					value(iter.actor());

				if (this->lpConstantCovariate)
				{
					alterValue *=
						this->lpConstantCovariate->value(iter.actor());
				}
				else if (this->lpChangingCovariate)
				{
					alterValue *=
						this->lpChangingCovariate->value(iter.actor(),
							period);
				}
				else
					throw logic_error(
						"No individual covariate found.");

				totalAlterValue += alterValue;
			}
		}
		if (fabs(totalAlterValue) < 1e-6)
		{
			return 1;
		}
		else
		{
			return pow(this->lpTable->value(1,1), totalAlterValue);
// seems strange; can't this be made more efficient?
		}
	}
	else
	{
		throw new logic_error("Unexpected diffusion rate effect type"+this->leffectName);
	}
}
/**
 * Stores the parameter for the diffusion rate effect.
 */
void DiffusionRateEffect::parameter(double parameterValue) const
{
	this->lpTable->parameter(parameterValue);
}

/**
 * Returns the parameter for the diffusion rate effect.
 */
double DiffusionRateEffect::parameter() const
{
	return this->lpTable->parameter();
}

/**
 * Stores the internal effect parameter for the diffusion rate effect.
 */
void DiffusionRateEffect::internalEffectParameter(int parValue)
{
	this->linternalEffectParameter = parValue;
	this->linternalNonZero = (this->linternalEffectParameter != 0);
}

/**
 * Returns the internal effect parameter for the diffusion rate effect.
 */
int DiffusionRateEffect::internalEffectParameter() const
{
	return this->linternalEffectParameter;
}


}
