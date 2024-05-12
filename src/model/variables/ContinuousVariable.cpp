/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ContinuousVariable.cpp
 *
 * Description: This file contains the implementation of the
 * ContinuousVariable class.
 *****************************************************************************/

#include <R_ext/Arith.h>
#include <R_ext/Print.h>
#include <R_ext/Error.h>

#include "ContinuousVariable.h"
#include "data/ContinuousLongitudinalData.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "model/Model.h"
#include "model/SdeSimulation.h"
#include "model/SimulationActorSet.h"
#include "model/effects/ContinuousEffect.h"
#include "model/effects/EffectFactory.h"

using namespace std;

namespace siena
{
/**
 * Creates a new continuous behavior variable for the given observed data.
 * @param pSimulation the model simulation, which becomes the owner of
 * this variable
 */
ContinuousVariable::ContinuousVariable(ContinuousLongitudinalData * pData,
	EpochSimulation * pSimulation) : NamedObject(pData->name())
{
	this->lpActorSet = pSimulation->pSimulationActorSet(pData->pActorSet());
	this->lpSimulation = pSimulation;

	this->lpData = pData;
	this->lvalues = new double[this->n()];
	this->lpFunction = new Function();

	this->leffectContribution = new double * [this->n()];
	for (int i = 0; i < this->n(); i++)
	{
		this->leffectContribution[i] =
			new double[pSimulation->pModel()->rEvaluationEffects(pData->name()).size()];
	}
}

/**
 * Deallocates this variable object.
 */
ContinuousVariable::~ContinuousVariable()
{
	for (int i = 0; i < this->n(); i++)
	{
		delete[] this->leffectContribution[i];
	}
	delete[] this->leffectContribution;
	delete this->lpFunction;
	delete[] this->lvalues;

	this->lpActorSet = 0;
	this->lpSimulation = 0;
	this->lpData = 0;
	this->leffectContribution = 0;
	this->lvalues = 0;
	this->lpFunction = 0;
}


// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Initializes this dependent variable as of the beginning of the given period.
 */
void ContinuousVariable::initialize(int period)
{
	this->lperiod = period;
	this->lsimulatedDistance = 0;

	this->lbasicScale = this->lpSimulation->pModel()->basicScaleParameter(period);
	this->lbasicScaleScore = 0;
	this->lbasicScaleDerivative = 0;

	// Copy the values from the corresponding observation.
	for (int i = 0; i < this->n(); i++)
	{
		this->lvalues[i] = this->lpData->value(period, i);
	}
}

/**
 * Creates SDE effects and stores them in a function.
 */
void ContinuousVariable::initializeFunction() const
{
	const vector<EffectInfo *> & rEffects =
		this->lpSimulation->pModel()->rEvaluationEffects(this->name());

	EffectFactory factory(this->lpSimulation->pData());

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pEffectInfo = rEffects[i];
		Effect * pEffect = factory.createEffect(pEffectInfo);
		this->lpFunction->addEffect(pEffect);
	}
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the number of actors for this variable.
 */
int ContinuousVariable::n() const
{
	return this->lpActorSet->n();
}


/**
 * Returns the ID of this dependent variable, which is the same as the ID
 * of the underlying observed data object.
 */
int ContinuousVariable::id() const
{
	return this->lpData->id();
}


/**
 * Stores the distance of this variable to the observed data at the
 * beginning of the current period.
 */
void ContinuousVariable::simulatedDistance(double distance)
{
	this->lsimulatedDistance = distance;
}


/**
 * Returns the distance of this variable to the observed data at the
 * beginning of the current period.
 */
double ContinuousVariable::simulatedDistance() const
{
	return this->lsimulatedDistance;
}


/**
 * Returns the current value on this variable for the given actor.
 */
double ContinuousVariable::value(int actor) const
{
	return this->lvalues[actor];
}


/**
 * Stores the current value on this variable for the given actor.
 */
void ContinuousVariable::value(int actor, double newValue)
{
	this->lvalues[actor] = newValue;
}


// ----------------------------------------------------------------------------
// Section: Changing the behavior variable
// ----------------------------------------------------------------------------


/**
 * Calculates the contribution of all individual effects per actor and
 * stores them in a 2D array.
 */
void ContinuousVariable::calculateEffectContribution()
{
	const Function * pFunction = this->pFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		ContinuousEffect * pEffect =
			(ContinuousEffect *) pFunction->rEffects()[i];

		for (int actor = 0; actor < this->n(); actor++)
		{
			this->leffectContribution[actor][i] =
				pEffect->calculateChangeContribution(actor);
		}
	}
}

/**
 * Returns the total contribution of all effects in the SDE for a certain
 * actor, except for the random effects.
 */
double ContinuousVariable::totalFunctionContribution(int actor) const
{
	double contribution = 0;
	const Function * pFunction = this->pFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		ContinuousEffect * pEffect =
			(ContinuousEffect *) pFunction->rEffects()[i];
        contribution += pEffect->coefficient() * leffectContribution[actor][i];
	}

	return contribution;
}


/**
 * Updates the scores for effects according to the current step in the
 * simulation.
 */
void ContinuousVariable::accumulateScores(const vector<double> &actorMeans,
	                                const vector<double> &actorErrors,
									double dt)
{
	const SdeSimulation * sde = this->pSimulation()->pSdeSimulation();
	double a = sde->feedbackParameter();
	double Adt = sde->feedbackCoefficient();
	double g = sde->wienerParameter();
	double Vardt = sde->wienerCoefficient();
	double tau = this->pSimulation()->pModel()->basicScaleParameter(
										this->pSimulation()->period());

	double errorxerror = 0; // inner product
	for (int actor = 0; actor < this->n(); actor++)
	{
		errorxerror += actorErrors[actor] * actorErrors[actor];
	}

	vector<double> bxeffects(this->n()); // inner product
	for (int actor = 0; actor < this->n(); actor++) 
	{
		bxeffects[actor]= 0;
	}
	for (unsigned i = 0; i < this->pFunction()->rEffects().size(); i++)
	{
		Effect * pEffect = this->pFunction()->rEffects()[i];
		if (pEffect->pEffectInfo()->effectName() != "feedback" &&
			pEffect->pEffectInfo()->effectName() != "wiener")
		{
			for (int actor = 0; actor < this->n(); actor++)
			{
				bxeffects[actor] += pEffect->parameter() * leffectContribution[actor][i];
			}
		}
	}

	for (unsigned i = 0; i < this->pFunction()->rEffects().size(); i++)
	{
		Effect * pEffect = this->pFunction()->rEffects()[i];
		double score;

		if (pEffect->pEffectInfo()->effectName() == "feedback")
		{
			score = this->n() / (2*a) * (1 - g*g*tau*dt*Adt*Adt / Vardt);
			double C1 = 1 / (2*a*Vardt) * (1 - tau*dt*g*g*Adt*Adt / Vardt);
			double C2 = 0;
			for (int actor = 0; actor < this->n(); actor++)
			{
				double temp = tau*dt*actorMeans[actor] + bxeffects[actor] / a * (tau*dt - (Adt-1)/a);
				C2 += actorErrors[actor] * temp;
			}
			C2 *= -2;

			score += - C1 * errorxerror - 1/(2*Vardt) * C2;
		}
		else if (pEffect->pEffectInfo()->effectName() == "wiener")
		{
			score = -this->n() / g  + 1 / (g * Vardt) * errorxerror;
		}
		else // scores for parameters b
		{
			double errorxeffect = 0; // inner product
			for (int actor = 0; actor < this->n(); actor++)
			{
				errorxeffect += actorErrors[actor] * leffectContribution[actor][i];
			}
			score = 2 / ((Adt + 1)*g*g) * errorxeffect;
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}

	// score for parameter tau
	double score = - this->n() * g*g*dt*Adt*Adt / (2*Vardt);
	double C1 = - g*g*dt*Adt*Adt / (2*Vardt*Vardt);
	double C2 = 0;
	for (int actor = 0; actor < this->n(); actor++)
	{
		double temp = a*actorMeans[actor] + bxeffects[actor];
		C2 += actorErrors[actor] * temp;
	}
	C2 *= -2*dt;

	score += - C1 * errorxerror - 1/(2*Vardt) * C2;

	score += this->pSimulation()->pSdeSimulation()->basicScaleScore();
	this->pSimulation()->basicScaleScore(score);
}

}
