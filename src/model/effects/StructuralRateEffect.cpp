/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: StructuralRateEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * StructuralRateEffect.
 *****************************************************************************/

#include "StructuralRateEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/EffectValueTable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pVariable the network variable this effect depends on
 * @param[in] type the type of this effect
 * @param[in] parameter the statistical parameter of this effect
 */
StructuralRateEffect::StructuralRateEffect(const NetworkVariable * pVariable,
	StructuralRateEffectType type,
	double parameter)
{
	this->lpVariable = pVariable;
	this->ltype = type;

	double possibleDegree = std::max(this->lpVariable->n(),
		this->lpVariable->m());
	if (this->ltype == INVERSE_OUT_DEGREE_RATE)
	{
		this->lpTable = new EffectValueTable(possibleDegree, invertor);
	}
	else if (this->ltype == LOG_OUT_DEGREE_RATE)
	{
		this->lpTable = new EffectValueTable(possibleDegree, logarithmer);
	}
	else
	{
		this->lpTable = new EffectValueTable(possibleDegree, identity);
	}

	this->lpTable->parameter(parameter);
}

/**
 * Destructor.
 */
StructuralRateEffect::~StructuralRateEffect()
{
	delete this->lpTable;
	this->lpTable = 0;
}

/**
 * Returns the contribution of this effect for the given actor.
 */
double StructuralRateEffect::value(int i) const
{
	Network * pNetwork = this->lpVariable->pNetwork();

	switch (this->ltype)
	{
		case OUT_DEGREE_RATE:
		case INVERSE_OUT_DEGREE_RATE:
		case LOG_OUT_DEGREE_RATE:
			return this->lpTable->value(pNetwork->outDegree(i));

		case IN_DEGREE_RATE:
			return this->lpTable->value(pNetwork->inDegree(i));

		case RECIPROCAL_DEGREE_RATE:
			return this->lpTable->value(
				(((OneModeNetwork *) pNetwork)->reciprocalDegree(i)));
	}

	throw std::logic_error("Unexpected structural rate effect type");
}
/**
 * Stores the parameter for the structural rate effect.
 */
void StructuralRateEffect::parameter(double parameterValue)
{
	this->lpTable->parameter(parameterValue);
}

/**
 * Returns the parameter for the structural rate effect.
 */
double StructuralRateEffect::parameter() const
{
	return this->lpTable->parameter();
}

}
