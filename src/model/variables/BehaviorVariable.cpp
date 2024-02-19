/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorVariable.cpp
 *
 * Description: This file contains the implementation of the
 * BehaviorVariable class.
 *****************************************************************************/
#include <cstdlib>
#include <cmath>
#include <string>
#include <stdexcept>
#include <R_ext/Print.h>
#include <R_ext/Arith.h>
#include <R_ext/Error.h>
#include "data/ActorSet.h"
#include "utils/Random.h"
#include <Rinternals.h>
#include <R_ext/Print.h>
#include <R_ext/Arith.h>
#include "BehaviorVariable.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/tables/Cache.h"
#include "DependentVariable.h"
#include "model/Model.h"
#include "model/effects/BehaviorEffect.h"
#include "model/EffectInfo.h"
#include "model/SimulationActorSet.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/ml/BehaviorChange.h"
#include <Rinternals.h>

using namespace std;

namespace siena
{
SEXP getMiniStepDF(const MiniStep& miniStep);

/**
 * Creates a new behavior variable for the given observed data.
 * @param pModel the owner model of this variable
 */
BehaviorVariable::BehaviorVariable(BehaviorLongitudinalData * pData,
	EpochSimulation * pSimulation) :
		DependentVariable(pData->name(),
			pData->pActorSet(),
			pSimulation)
{
	this->lpData = pData;
	this->lvalues = new int[this->n()];
	this->levaluationEffectContribution = new double * [3];
	this->lendowmentEffectContribution = new double * [3];
	this->lcreationEffectContribution = new double * [3];
	this->lprobabilities = new double[3];

	for (int i = 0; i < 3; i++)
	{
		this->levaluationEffectContribution[i] =
			new double[pSimulation->pModel()->rEvaluationEffects(pData->name()).size()];
		this->lendowmentEffectContribution[i] =
			new double[pSimulation->pModel()->rEndowmentEffects(pData->name()).size()];
		this->lcreationEffectContribution[i] =
			new double[pSimulation->pModel()->rCreationEffects(pData->name()).size()];
		this->lprobabilities[i] = 0;
	}

	this->lbehaviorModelType = BehaviorModelType(pData->behModelType());
	this->lego = 0;
	this->lupPossible = true;
	this->ldownPossible = true;
}


/**
 * Deallocates this variable object.
 */
BehaviorVariable::~BehaviorVariable()
{
	delete[] this->lvalues;

	this->lpData = 0;
	this->lvalues = 0;
	delete[] this->lprobabilities;
	// Delete arrays of contributions

	for (int i = 0; i < 3; i++)
	{
		delete[] this->levaluationEffectContribution[i];
		delete[] this->lendowmentEffectContribution[i];
		delete[] this->lcreationEffectContribution[i];
	}

	delete[] this->levaluationEffectContribution;
	delete[] this->lendowmentEffectContribution;
	delete[] this->lcreationEffectContribution;

	this->levaluationEffectContribution = 0;
	this->lendowmentEffectContribution = 0;
	this->lcreationEffectContribution = 0;
	this->lprobabilities = 0;

	// no need to delete lpChangeContribution since this is
	// handled by the MiniStep
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the second dimension of this variable, namely, how many values
 * correspond to each actor. This number is 1 for behavior variables.
 */
int BehaviorVariable::m() const
{
	return 1;
}


/**
 * Returns the longitudinal data object this variable is based on.
 */
LongitudinalData * BehaviorVariable::pData() const
{
	return this->lpData;
}


/**
 * Returns if the observed start value of the given actor is missing
 * for the current period.
 */
bool BehaviorVariable::missingStartValue(int actor) const
{
	return this->lpData->missing(this->period(), actor);
}


/**
 * Returns the current value on this behavior for the given actor.
 */
int BehaviorVariable::value(int actor) const
{
	return this->lvalues[actor];
}


/**
 * Stores the current value on this behavior for the given actor.
 */
void BehaviorVariable::value(int actor, int newValue)
{
	this->lvalues[actor] = newValue;
}


/**
 * Returns the current value on this behavior for the given actor, which is
 * centered around the overall mean of the observed values.
 */
double BehaviorVariable::centeredValue(int actor) const
{
	return this->lvalues[actor] - this->lpData->overallMean();
}


/**
 * Returns if the behavior is structurally determined for the given actor
 * at the current period.
 */
bool BehaviorVariable::structural(int actor) const
{
	return this->lpData->structural(this->period(), actor);
}


/**
 * Returns the centered similarity of the given actors.
 */
double BehaviorVariable::similarity(int i, int j) const
{
	return this->lpData->similarity(this->lvalues[i], this->lvalues[j]);
}


/**
 * Returns the array of current values for this behavior variable.
 */
const int * BehaviorVariable::values() const
{
	return this->lvalues;
}


/**
 * Returns the range of observed values for this behavior variable.
 */
int BehaviorVariable::range() const
{
	return this->lpData->range();
}

/**
 * Returns the similarity mean for this behavior variable.
 */
double BehaviorVariable::similarityMean() const
{
	return this->lpData->similarityMean();
}

/**
 * Sets the model type.
 */
void BehaviorVariable::behaviorModelType(int type)
{
	this->lbehaviorModelType = BehaviorModelType(type);
}

/**
 * Returns the model type.
 */
BehaviorModelType BehaviorVariable::behaviorModelType() const
{
	return this->lbehaviorModelType;
}

/**
 * Returns whether the model type is absorb.
 */

bool BehaviorVariable::behaviorModelTypeABSORB() const
{
	return (this->lbehaviorModelType == ABSORB);
}



// ----------------------------------------------------------------------------
// Section: Initialization at the beginning of a period
// ----------------------------------------------------------------------------

/**
 * Initializes this variable as of the beginning of the given period.
 */
void BehaviorVariable::initialize(int period)
{
	DependentVariable::initialize(period);

	// Copy the values from the corresponding observation.

	for (int i = 0; i < this->n(); i++)
	{
		this->lvalues[i] = this->lpData->value(period, i);
	}

	this->behaviorModelType(this->lpData->behModelType());
}


// ----------------------------------------------------------------------------
// Section: Composition change
// ----------------------------------------------------------------------------

/**
 * Sets leavers values back to the value at the start of the simulation.
 */
void BehaviorVariable::setLeaverBack(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->pActorSet())
	{
		// Reset ties from the given actor to values at start

		for (int i = 0; i < this->n(); i++)
		{
			this->lvalues[actor] =	this->lpData->value(this->period(), actor);
		}
	}
}


// ----------------------------------------------------------------------------
// Section: Changing the behavior variable
// ----------------------------------------------------------------------------

/**
 * Simulates a change of the behavior according to the choice of the given
 * actor.
 */
void BehaviorVariable::makeChange(int actor)
{
	this->lego = actor;
	this->calculateProbabilities(actor);

	// Choose the change
	int difference = nextIntWithProbabilities(3, this->lprobabilities) - 1;

	if (difference < -1)
	{
		difference = -1;
	}
	else if (difference > 1)
	{
		difference = 1;
	}

	if (this->pSimulation()->pModel()->needScores())
	{
		this->accumulateScores(difference + 1,
			this->lupPossible,
			this->ldownPossible);
	}
	if (this->pSimulation()->pModel()->needDerivatives())
	{
		this->accumulateDerivatives(); // ABC
	}

	if (this->pSimulation()->pModel()->needChain())
	{
		// insert ministep in chain
		BehaviorChange * pMiniStep =
			new BehaviorChange(this->lpData, actor, difference);
		if (this->pSimulation()->pModel()->needChangeContributions())
		{
			pMiniStep->changeContributions(lpChangeContribution);
		}
		this->pSimulation()->pChain()->insertBefore(pMiniStep,
			this->pSimulation()->pChain()->pLast());
		pMiniStep->logChoiceProbability(log(this->lprobabilities[difference
					+ 1]));
	}
	// Make the change

	if (difference != 0)
	{
		int oldValue = this->lvalues[actor];

		// Make the change
		this->lvalues[actor] += difference;

		// Update the distance from the observed data at the beginning of the
		// period. Actors with missing values at any of the endpoints of the
		// period don't contribute to the distance

		if (!this->lpData->missing(this->period(), actor) &&
			!this->lpData->missing(this->period() + 1, actor))
		{
			int observedValue = this->lpData->value(this->period(), actor);
			this->simulatedDistance(this->simulatedDistance() +
				std::abs(this->lvalues[actor] - observedValue) -
				std::abs(oldValue - observedValue));
		}
	}
	this->successfulChange(true);
}


/**
 * Calculates the probabilities of each possible change.
 */
void BehaviorVariable::calculateProbabilities(int actor)
{
	double maxValue = 0.0; // used to be R_NegInf, but always there is a 0

	this->preprocessEgo();
	this->lupPossible = true;
	this->ldownPossible = true;
	int currentValue = this->lvalues[actor];
	bool ismax = (currentValue >= this->lpData->max());
	bool ismin = (currentValue <= this->lpData->min());

	unsigned evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	unsigned endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();
	unsigned creationEffectCount = this->pCreationFunction()->rEffects().size();

	// initialize for later use!
	for (unsigned i = 0; i < evaluationEffectCount; i++)
	{
		this->levaluationEffectContribution[1][i] =	0;
	}
	for (unsigned i = 0; i < endowmentEffectCount; i++)
	{
		this->lendowmentEffectContribution[1][i] =	0;
		this->lendowmentEffectContribution[2][i] = 	0;
	}

	for (unsigned i = 0; i < creationEffectCount; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->lcreationEffectContribution[j][i] = 0;
		}
	}

	if (this->pSimulation()->pModel()->needChangeContributions())
	{
		this->lpChangeContribution = new map<const EffectInfo *, vector<double> >();
// used to be:   >[evaluationEffectCount+endowmentEffectCount+creationEffectCount];
		for (unsigned i = 0; i < evaluationEffectCount; i++)
		{
			vector<double> vec(3,0);
			this->lpChangeContribution->insert(make_pair(this->pEvaluationFunction()->rEffects()[i]->pEffectInfo(), vec));
		}
		for (unsigned i = 0; i < endowmentEffectCount; i++)
		{
			vector<double> vec(3,0);
			this->lpChangeContribution->insert(make_pair(this->pEndowmentFunction()->rEffects()[i]->pEffectInfo(), vec));
		}
		for (unsigned i = 0; i < creationEffectCount; i++)
		{
			vector<double> vec(3,0);
			this->lpChangeContribution->insert(make_pair(this->pCreationFunction()->rEffects()[i]->pEffectInfo(), vec));
		}
	}

	// Calculate the objective function for downward change.
	// Defer exp until we can subtract the largest to avoid overflow.
	if (ismin || this->lpData->upOnly(this->period()))
	{
		this->lprobabilities[0] = 0;
		this->ldownPossible = false;
		for (unsigned i = 0; i < pEvaluationFunction()->rEffects().size(); i++)
		{
			this->levaluationEffectContribution[0][i] =	R_NaN;
		}
		for (unsigned i = 0; i < pEndowmentFunction()->rEffects().size(); i++)
		{
			this->lendowmentEffectContribution[0][i] =  R_NaN;
		}

		for (unsigned i = 0;
			i < this->pCreationFunction()->rEffects().size();
			i++)
		{
			this->lcreationEffectContribution[0][i] = R_NaN;
		}
	}
	else
	{
		this->lprobabilities[0] =
			this->totalEvaluationContribution(actor, -1) +
				this->totalEndowmentContribution(actor, -1);
		maxValue = max(maxValue, this->lprobabilities[0]);
	}

	// No change means zero contribution, but exp(0) = 1
	this->lprobabilities[1] = 0;

	// Calculate the objective function for upward change
	if (ismax || this->lpData->downOnly(this->period()))
	{
		this->lprobabilities[2] = 0;
		this->lupPossible = false;
		for (unsigned i = 0; i < pEvaluationFunction()->rEffects().size(); i++)
		{
			this->levaluationEffectContribution[2][i] =	R_NaN;
		}
		for (unsigned i = 0; i < pEndowmentFunction()->rEffects().size(); i++)
		{
			this->lendowmentEffectContribution[2][i] = R_NaN;
		}
		for (unsigned i = 0;
				i < this->pCreationFunction()->rEffects().size();
				i++)
		{
			this->lcreationEffectContribution[2][i] = R_NaN;
		}
	}
	else
	{
		this->lprobabilities[2] =
			this->totalEvaluationContribution(actor, 1) +
				this->totalCreationContribution(actor, 1);
		maxValue = max(maxValue, this->lprobabilities[2]);
	}

		// turn objective function values into probabilities
	if (this->ldownPossible)
	{
		this->lprobabilities[0] = exp(this->lprobabilities[0] - maxValue);
	}
	if (this->lupPossible)
	{
		this->lprobabilities[2] = exp(this->lprobabilities[2] - maxValue);
	}
	this->lprobabilities[1]  = exp(-maxValue);

	double sum = 0;
	if ((this->behaviorModelTypeABSORB()) && (ismin || ismax))
	{
		if (ismin)
		{
			sum = 2 * this->lprobabilities[1] + this->lprobabilities[2];
			this->lprobabilities[1] = 2 * this->lprobabilities[1]/sum;
			this->lprobabilities[2] /= sum;
		}
		else
		{
			sum = 2 * this->lprobabilities[1] + this->lprobabilities[0];
			this->lprobabilities[1] = 2 * this->lprobabilities[1]/sum;
			this->lprobabilities[0] /= sum;
		}
	}
	else
	{
		sum = this->lprobabilities[0] + this->lprobabilities[1] +
							this->lprobabilities[2];
		this->lprobabilities[0] /= sum;
		this->lprobabilities[1] /= sum;
		this->lprobabilities[2] /= sum;
	}
}


/**
 * Returns the total contribution of all effects in the evaluation function if
 * the behavior of the given actor is changed by the given amount.
 */
double BehaviorVariable::totalEvaluationContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		BehaviorEffect * pEffect = (BehaviorEffect *) pFunction->rEffects()[i];
		double thisContribution =
			pEffect->calculateChangeContribution(actor, difference);
		if (this->pSimulation()->pModel()->needChangeContributions())
		{
			(* this->lpChangeContribution)[pEffect->pEffectInfo()]
				.at(difference + 1) = thisContribution;
		}
		this->levaluationEffectContribution[difference + 1][i] =
			thisContribution;
		contribution += pEffect->parameter() * thisContribution;
	}
	return contribution;
}

double BehaviorVariable::totalEndowmentContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Function * pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		BehaviorEffect * pEffect = (BehaviorEffect *) pFunction->rEffects()[i];
		double thisContribution =
			pEffect->calculateChangeContribution(actor, difference);
		if (this->pSimulation()->pModel()->needChangeContributions())
		{
			(* this->lpChangeContribution)[pEffect->pEffectInfo()]
				.at(difference + 1) = thisContribution;
		}
		this->lendowmentEffectContribution[difference + 1][i] =
			thisContribution;
		contribution += pEffect->parameter() * thisContribution;
	}

	return contribution;
}


/**
 * Returns the total contribution of all effects in the tie creation function
 * if the behavior of the given actor is changed by the given amount.
 */
double BehaviorVariable::totalCreationContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Function * pFunction = this->pCreationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		BehaviorEffect * pEffect = (BehaviorEffect *) pFunction->rEffects()[i];
		double thisContribution =
			pEffect->calculateChangeContribution(actor, difference);
		if (this->pSimulation()->pModel()->needChangeContributions())
		{
			(* this->lpChangeContribution)[pEffect->pEffectInfo()]
				.at(difference + 1) = thisContribution;
		}
		this->lcreationEffectContribution[difference + 1][i] = thisContribution;
		contribution += pEffect->parameter() * thisContribution;
	}

	return contribution;
}


/**
 * Updates the scores for effects according
 * to the current step in the simulation.
 * This function is called with arguments difference = 0, 1, 2.
 * 1 means no change.
 */
void BehaviorVariable::accumulateScores(int difference,
	bool upPossible, bool downPossible) const
{
	for (unsigned i = 0;
		i < this->pEvaluationFunction()->rEffects().size();
		i++)
	{
// 		if (difference == 1) no change, but not initialised
// 		{
// 			this->levaluationEffectContribution[difference][i] = 0;
// 		}
		Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
		double score = this->levaluationEffectContribution[difference][i];
		if (upPossible)
		{
			score -=
				this->levaluationEffectContribution[2][i] *
				this->lprobabilities[2];
		}

		if (downPossible)
		{
			score -=
				this->levaluationEffectContribution[0][i] *
				this->lprobabilities[0];
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
		if (R_IsNaN(score))
		{
			Rf_error("nan in accumulateScores1");
		}
	}

	for (unsigned i = 0;
		i < this->pEndowmentFunction()->rEffects().size();
		i++)
	{
// 		if (difference == 1) // no change, but not initialised
// 		{
// 			this->lendowmentEffectContribution[difference][i] = 0;
// 		}
// 		if (difference == 2) // up has no effect on endowment function
// 		{
// 			this->lendowmentEffectContribution[difference][i] = 0;
// 		}
		Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];
		double score = 0;

		// if difference is 1, no change so contribution is 0.
		if (difference == 0)
		{
			score = this->lendowmentEffectContribution[difference][i];
		}

		if (downPossible)
		{
			score -=
				this->lendowmentEffectContribution[0][i] *
					this->lprobabilities[0];

		}
		if (R_IsNaN(score))
		{
			Rf_error("nan in accumulateScores2");
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}

	for (unsigned i = 0;
		i < this->pCreationFunction()->rEffects().size();
		i++)
	{
		Effect * pEffect = this->pCreationFunction()->rEffects()[i];
		double score = 0;

		// if difference is 1, no change so contribution is 0.
		if (difference == 2)
		{
			score = this->lcreationEffectContribution[difference][i];
		}

		if (upPossible)
		{
			score -=
				this->lcreationEffectContribution[2][i] *
					this->lprobabilities[2];
		}

		if (R_IsNaN(score))
		{
			Rf_error("nan in accumulateScores3");
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}
}

/**
 * This method does some preprocessing to enable subsequent queries regarding
 * the current ego.
 */
void BehaviorVariable::preprocessEgo()
{
	// Let the effects do their preprocessing.

	this->preprocessEffects(this->pEvaluationFunction());
	this->preprocessEffects(this->pEndowmentFunction());
	this->preprocessEffects(this->pCreationFunction());
}


/**
 * This method does some preprocessing for each effect in the given function
 * to enable subsequent queries regarding the current ego.
 */
void BehaviorVariable::preprocessEffects(const Function * pFunction)
{
	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		BehaviorEffect * pEffect =
			(BehaviorEffect *) pFunction->rEffects()[i];
		pEffect->preprocessEgo(this->lego);
	}
}


// ----------------------------------------------------------------------------
// Section: Maximum likelihood related methods
// ----------------------------------------------------------------------------

/**
 * Calculates the probability of the given ministep assuming that the
 * ego of the ministep will change this variable.
 */
double BehaviorVariable::probability(MiniStep * pMiniStep)
{
	// Initialize the cache object for the current ego
	this->pSimulation()->pCache()->initialize(pMiniStep->ego());

	BehaviorChange * pBehaviorChange =
		dynamic_cast<BehaviorChange *>(pMiniStep);

	if (pBehaviorChange->difference() < -1 ||
		pBehaviorChange->difference() > 1)
	{
		throw invalid_argument("MiniStep difference out of range [-1,1].");
	}

	this->calculateProbabilities(pMiniStep->ego());
	if (this->pSimulation()->pModel()->needScores())
	{
		this->accumulateScores(pBehaviorChange->difference() + 1,
				this->lupPossible, this->ldownPossible);
	}
	if (this->pSimulation()->pModel()->needDerivatives())
	{
		this->accumulateDerivatives();
	}
	return this->lprobabilities[pBehaviorChange->difference() + 1];
}

/**
 * Updates the derivatives for effects
 * according to the current miniStep in the chain.
 */
void BehaviorVariable::accumulateDerivatives() const
{
	int totalEvaluationEffects = this->pEvaluationFunction()->rEffects().size();
	int totalEndowmentEffects = this->pEndowmentFunction()->rEffects().size();
	int totalCreationEffects = this->pCreationFunction()->rEffects().size();
	int totalEffects =
		totalEvaluationEffects + totalEndowmentEffects + totalCreationEffects;
	int firstCreationIndex = totalEvaluationEffects + totalEndowmentEffects;
	Effect * pEffect1;
	Effect * pEffect2;
	double derivative;
	double * product = new double[totalEffects];
	double contribution1 = 0.0;
	double contribution2 = 0.0;

	for (int effect1 = 0; effect1 < totalEffects; effect1++)
	{
		product[effect1] = 0.0;
		int endowment1 = effect1 - totalEvaluationEffects;
		int creation1 = effect1 - firstCreationIndex;

		if (effect1 < totalEvaluationEffects)
		{
			pEffect1 = this->pEvaluationFunction()->rEffects()[effect1];
		}
		else if (effect1 < firstCreationIndex)
		{
			pEffect1 = this->pEndowmentFunction()->rEffects()[endowment1];
		}
		else
		{
			pEffect1 = this->pCreationFunction()->rEffects()[creation1];
		}

		if (this->lupPossible)
		{
			if (effect1 < totalEvaluationEffects)
			{
				product[effect1] +=
					this->levaluationEffectContribution[2][effect1] *
					this->lprobabilities[2];
				if (R_IsNaN(product[effect1]))
				{
					Rprintf("eval up effect 1 %d %f %f \n", effect1,
						this->levaluationEffectContribution[2][effect1],
						this->lprobabilities[2]);
				}
			}
			else if (effect1 >= firstCreationIndex)
			{
				product[effect1] +=
					this->lcreationEffectContribution[2][creation1] *
					this->lprobabilities[2];

				if (R_IsNaN(product[effect1]))
				{
					Rprintf("creation up effect 1 %d %f %f \n",
						effect1,
						this->lcreationEffectContribution[2][creation1],
						this->lprobabilities[2]);
				}
			}
		}
		if (this->ldownPossible)
		{
			if (effect1 < totalEvaluationEffects)
			{
				product[effect1] +=
					this->levaluationEffectContribution[0][effect1] *
					this->lprobabilities[0];
				if (R_IsNaN(product[effect1]))
				{
					Rprintf("eval down effect 1 %d %f %f \n", effect1,
						this->levaluationEffectContribution[0][effect1],
						this->lprobabilities[0]);
				}
			}
			else if (effect1 < firstCreationIndex)
			{
				product[effect1] +=
					this->lendowmentEffectContribution[0][endowment1] *
					this->lprobabilities[0];
				if (R_IsNaN(product[effect1]))
				{
					Rprintf("endow down effect 1 %d %d%f %f \n", effect1,
						endowment1,
						this->lendowmentEffectContribution[0][endowment1],
						this->lprobabilities[0]);
				}
			}
		}
		for (int effect2 = effect1; effect2 < totalEffects; effect2++)
		{
			int endowment2 = effect2 - totalEvaluationEffects;
			int creation2 = effect2 - firstCreationIndex;

			derivative = 0.0;

			if (effect2 < totalEvaluationEffects)
			{
				pEffect2 = this->pEvaluationFunction()->rEffects()[effect2];
			}
			else if (effect2 < firstCreationIndex)
			{
				pEffect2 = this->pEndowmentFunction()->rEffects()[endowment2];
			}
			else
			{
				pEffect2 = this->pCreationFunction()->rEffects()[creation2];
			}

			if (this->ldownPossible &&
				effect1 < firstCreationIndex &&
				effect2 < firstCreationIndex)
			{
				if (effect1 < totalEvaluationEffects)
				{
					contribution1 =
						this->levaluationEffectContribution[0][effect1];
				}
				else
				{
					contribution1 =
						this->lendowmentEffectContribution[0][endowment1];
				}

				if (effect2 < totalEvaluationEffects)
				{
					contribution2 =
						this->levaluationEffectContribution[0][effect2];
				}
				else
				{
					contribution2 =
						this->lendowmentEffectContribution[0][endowment2];
				}

				derivative -=
					contribution1 * contribution2 *	this->lprobabilities[0];
			}

			if (this->lupPossible &&
				(effect1 < totalEvaluationEffects ||
					effect1 >= firstCreationIndex) &&
				(effect2 < totalEvaluationEffects ||
					effect2 >= firstCreationIndex))
			{
				if (effect1 < totalEvaluationEffects)
				{
					contribution1 =
						this->levaluationEffectContribution[2][effect1];
				}
				else
				{
					contribution1 =
						this->lcreationEffectContribution[2][creation1];
				}

				if (effect2 < totalEvaluationEffects)
				{
					contribution2 =
						this->levaluationEffectContribution[2][effect2];
				}
				else
				{
					contribution2 =
						this->lcreationEffectContribution[2][creation2];
				}

				derivative -=
					contribution1 * contribution2 *	this->lprobabilities[2];
			}

			this->pSimulation()->derivative(pEffect1->pEffectInfo(),
				pEffect2->pEffectInfo(),
				this->pSimulation()->derivative(pEffect1->pEffectInfo(),
					pEffect2->pEffectInfo()) +	derivative);
		}
	}

	for (int effect1 = 0; effect1 < totalEffects; effect1++)
	{
		int endowment1 = effect1 - totalEvaluationEffects;
		int creation1 = effect1 - firstCreationIndex;

		for (int effect2 = effect1; effect2 < totalEffects; effect2++)
		{
			int endowment2 = effect2 - totalEvaluationEffects;
			int creation2 = effect2 - firstCreationIndex;

			if (effect1 < totalEvaluationEffects)
			{
				pEffect1 = this->pEvaluationFunction()->rEffects()[effect1];
			}
			else if (effect1 < firstCreationIndex)
			{
				pEffect1 = this->pEndowmentFunction()->rEffects()[endowment1];
			}
			else
			{
				pEffect1 = this->pCreationFunction()->rEffects()[creation1];
			}

			if (effect2 < totalEvaluationEffects)
			{
				pEffect2 = this->pEvaluationFunction()->rEffects()[effect2];
			}
			else if (effect2 < firstCreationIndex)
			{
				pEffect2 = this->pEndowmentFunction()->rEffects()[endowment2];
			}
			else
			{
				pEffect2 = this->pCreationFunction()->rEffects()[creation2];
			}

			if (R_IsNaN(product[effect1]))
			{
				Rprintf("effect 1 %d \n", effect1);
			}
			if (R_IsNaN(product[effect2]))
			{
				Rprintf("effect2 %d \n", effect2);
			}
			this->pSimulation()->derivative(pEffect1->pEffectInfo(),
				pEffect2->pEffectInfo(),
				this->pSimulation()->derivative(pEffect1->pEffectInfo(),
					pEffect2->pEffectInfo()) +
				product[effect1] * product[effect2]);
		}
	}
	delete[] product;
}

/**
 * Returns whether applying the given ministep on the current state of this
 * variable would be valid with respect to all constraints. One can disable
 * the checking of up-only and down-only conditions.
 */
bool BehaviorVariable::validMiniStep(const MiniStep * pMiniStep,
	bool checkUpOnlyDownOnlyConditions) const
{
	bool valid = DependentVariable::validMiniStep(pMiniStep);

	if (valid && !pMiniStep->diagonal())
	{
		const BehaviorChange * pBehaviorChange =
			dynamic_cast<const BehaviorChange *>(pMiniStep);
		int i = pMiniStep->ego();
		int d = pBehaviorChange->difference();
		int newValue = this->lvalues[i] + d;

		if (newValue < this->lpData->min() || newValue > this->lpData->max())
		{
			valid = false;
		}
		else if (checkUpOnlyDownOnlyConditions &&
			d > 0 &&
			this->lpData->downOnly(this->period()))
		{
			valid = false;
		}
		else if (checkUpOnlyDownOnlyConditions &&
			d < 0 &&
			this->lpData->upOnly(this->period()))
		{
			valid = false;
		}
		else
		{
			valid = !this->lpData->structural(this->period(), i);
		}
	}

	return valid;
}


/**
 * Generates a random ministep for the given ego.
 */
MiniStep * BehaviorVariable::randomMiniStep(int ego)
{
	this->pSimulation()->pCache()->initialize(ego);
	this->calculateProbabilities(ego);
	int difference = nextIntWithProbabilities(3, this->lprobabilities) - 1;
	BehaviorChange * pMiniStep =
		new BehaviorChange(this->lpData, ego, difference);
	pMiniStep->logChoiceProbability(log(this->lprobabilities[difference + 1]));
	return pMiniStep;
}


/**
 * Returns if the observed value for the option of the given ministep
 * is missing at either end of the period.
 */
bool BehaviorVariable::missing(const MiniStep * pMiniStep) const
{
	return this->lpData->missing(this->period(), pMiniStep->ego()) ||
		this->lpData->missing(this->period() + 1, pMiniStep->ego());
}

/**
 * Returns if the given ministep is structurally determined in the period.
 */
bool BehaviorVariable::structural(const MiniStep * pMiniStep) const
{
	return this->lpData->structural(this->period(), pMiniStep->ego());
}



// ----------------------------------------------------------------------------
// Section: Properties
// ----------------------------------------------------------------------------

/**
 * Returns if this is a behavior variable.
 */
bool BehaviorVariable::behaviorVariable() const
{
	return true;
}

}
