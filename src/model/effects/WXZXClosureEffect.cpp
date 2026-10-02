/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WXZXClosureEffect.cpp
 *
 * Description: Implementation of the WXZ => X multilevel four-cycle closure.
 * This class is a direct descendant of class NetworkEffect.
 * Alternatively, a class DoubleDyadicCovariateNetworkEffect could be created,
 * with WXZXClosureEffect as a descendant. 
 * This would be an option if more effects of this type would be created.
 *****************************************************************************/

#include <stdexcept>
#include <string>

#include "WXZXClosureEffect.h"
#include "data/ActorSet.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/Data.h"
#include "data/DyadicCovariate.h"
#include "data/DyadicCovariateValueIterator.h"
#include "model/EffectInfo.h"
#include "network/IncidentTieIterator.h"
#include "network/Network.h"

using namespace std;

namespace siena
{

WXZXClosureEffect::WXZXClosureEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
	this->lsums = 0;
	this->lexcludeMissings = false;
}

WXZXClosureEffect::~WXZXClosureEffect()
{
	delete[] this->lsums;
	this->lsums = 0;
}

void WXZXClosureEffect::initialize(const Data * pData,
	State * pState, int period, Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);

	string wName = this->pEffectInfo()->interactionName1();
	string zName = this->pEffectInfo()->interactionName2();

	this->lpWConstant = pData->pConstantDyadicCovariate(wName);
	this->lpWChanging = pData->pChangingDyadicCovariate(wName);
	this->lpZConstant = pData->pConstantDyadicCovariate(zName);
	this->lpZChanging = pData->pChangingDyadicCovariate(zName);

	if (!this->lpWConstant && !this->lpWChanging)
	{
		throw logic_error("First-mode dyadic covariate '" + wName +
			"' expected for WXZX.");
	}
	if (!this->lpZConstant && !this->lpZChanging)
	{
		throw logic_error("Second-mode dyadic covariate '" + zName +
			"' expected for WXZX.");
	}

	const DyadicCovariate * pW = this->lpWConstant ?
		static_cast<const DyadicCovariate *>(this->lpWConstant) :
		static_cast<const DyadicCovariate *>(this->lpWChanging);
	const DyadicCovariate * pZ = this->lpZConstant ?
		static_cast<const DyadicCovariate *>(this->lpZConstant) :
		static_cast<const DyadicCovariate *>(this->lpZChanging);

	int n = this->pNetwork()->n();
	int m = this->pNetwork()->m();

	if (pW->pFirstActorSet()->n() != n ||
		pW->pSecondActorSet()->n() != n)
	{
		throw logic_error("WXZX covariate W must be one-mode on the first "
			"node set of X.");
	}
	if (pZ->pFirstActorSet()->n() != m ||
		pZ->pSecondActorSet()->n() != m)
	{
		throw logic_error("WXZX covariate Z must be one-mode on the second "
			"node set of X.");
	}

	delete[] this->lsums;
	this->lsums = new double[m];
}

DyadicCovariateValueIterator WXZXClosureEffect::wRowValues(int i) const
{
	if (this->lpWConstant)
	{
		return this->lpWConstant->rowValues(i);
	}
	return this->lpWChanging->rowValues(i, this->period(),
		this->lexcludeMissings);
}

DyadicCovariateValueIterator WXZXClosureEffect::zRowValues(int k) const
{
	if (this->lpZConstant)
	{
		return this->lpZConstant->rowValues(k);
	}
	return this->lpZChanging->rowValues(k, this->period(),
		this->lexcludeMissings);
}

void WXZXClosureEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);

	int m = this->pNetwork()->m();
	for (int j = 0; j < m; j++)
	{
		this->lsums[j] = 0;
	}

	// i -W-> h -X-> k -Z-> j
	for (DyadicCovariateValueIterator iterH = this->wRowValues(ego);
		iterH.valid(); iterH.next())
	{
		int h = iterH.actor();
		double wih = iterH.value();

		for (IncidentTieIterator iterK = this->pNetwork()->outTies(h);
			iterK.valid(); iterK.next())
		{
			int k = iterK.actor();

			for (DyadicCovariateValueIterator iterJ = this->zRowValues(k);
				iterJ.valid(); iterJ.next())
			{
				int j = iterJ.actor();
				this->lsums[j] += wih * iterJ.value();
			}
		}
	}
}

double WXZXClosureEffect::calculateContribution(int alter) const
{
	return this->lsums[alter];
}

double WXZXClosureEffect::tieStatistic(int alter)
{
	return this->lsums[alter];
}

void WXZXClosureEffect::initializeStatisticCalculation()
{
	this->lexcludeMissings = true;
}

void WXZXClosureEffect::cleanupStatisticCalculation()
{
	this->lexcludeMissings = false;
}

}
