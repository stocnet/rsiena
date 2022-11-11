/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EffectFactory.cpp
 *
 * Description: This file contains the implementation of the
 * EffectFactory class.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>

#include "EffectFactory.h"
#include "data/Data.h"
#include "data/NetworkLongitudinalData.h"
#include "data/ContinuousLongitudinalData.h"
#include "model/EffectInfo.h"
#include "model/effects/AllEffects.h"
#include "model/effects/generic/GenericNetworkEffect.h"
#include "model/effects/generic/OutTieFunction.h"
#include "model/effects/generic/InTieFunction.h"
#include "model/effects/generic/ProductFunction.h"
#include "model/effects/generic/ConstantFunction.h"
#include "model/effects/generic/InDegreeFunction.h"
#include "model/effects/generic/IntSqrtFunction.h"
#include "model/effects/generic/DifferenceFunction.h"
#include "model/effects/generic/AbsDiffFunction.h"
#include "model/effects/generic/SumFunction.h"
#include "model/effects/generic/EgoFunction.h"
#include "model/effects/generic/EgoInDegreeFunction.h"
#include "model/effects/generic/OutDegreeFunction.h"
#include "model/effects/generic/EgoOutDegreeFunction.h"
#include "model/effects/generic/EgoRecipDegreeFunction.h"
#include "model/effects/generic/DegreeFunction.h"
#include "model/effects/generic/BetweennessFunction.h"
#include "model/effects/generic/GwespFunction.h"
#include "model/effects/generic/InStarFunction.h"
#include "model/effects/generic/OutStarFunction.h"
#include "model/effects/generic/InJaccardFunction.h"
#include "model/effects/generic/OutJaccardFunction.h"
#include "model/effects/generic/ReciprocatedTwoPathFunction.h"
#include "model/effects/generic/TwoPathFunction.h"
#include "model/effects/generic/TwoStepFunction.h"
#include "model/effects/generic/ReverseTwoPathFunction.h"
#include "model/effects/generic/MixedTwoStepFunction.h"
#include "model/effects/generic/WeightedMixedTwoPathFunction.h"
#include "model/effects/generic/ConditionalFunction.h"
#include "model/effects/generic/EqualCovariatePredicate.h"
#include "model/effects/generic/MissingCovariatePredicate.h"
#include "model/effects/generic/DoubleOutActFunction.h"
#include "model/effects/generic/CovariateDistance2AlterNetworkFunction.h"
#include "model/effects/generic/CovariateDistance2InAlterNetworkFunction.h"
#include "model/effects/generic/CovariateDistance2SimilarityNetworkFunction.h"
#include "model/effects/generic/CovariateDistance2EgoAltSimNetworkFunction.h"
#include "model/effects/generic/CovariateMixedNetworkAlterFunction.h"
#include "model/effects/generic/CovariateDegreeFunction.h"
#include "model/effects/generic/SameCovariateInStarFunction.h"
#include "model/effects/generic/SameCovariateOutStarFunction.h"
#include "model/effects/generic/SameCovariateTwoPathFunction.h"
#include "model/effects/generic/SameCovariateMixedTwoPathFunction.h"
#include "model/effects/generic/SameCovariateInTiesFunction.h"
#include "model/effects/generic/SameCovariateOutTiesFunction.h"
#include "model/effects/generic/DifferentCovariateInStarFunction.h"
#include "model/effects/generic/DifferentCovariateOutStarFunction.h"
#include "model/effects/generic/HomCovariateMixedTwoPathFunction.h"
#include "model/effects/generic/OutActDistance2Function.h"
#include "model/effects/generic/MixedThreeCyclesFunction.h"
#include "model/effects/generic/InStarsTimesDegreesFunction.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include "model/tables/NetworkCache.h"


using namespace std;

namespace siena
{
/**
 * GMM group information. Maps effect short names to the one unique gmm group
 * (usually the basic effect). In addition there are two special groups. The
 * empty string ("") meaning a match by position (diagonal matrix) and the
 * star ("*") meaning matching all.
 */
const std::string EffectFactory::gmmGroup(const EffectInfo* pEffectInfo) {
	// map basic rates to the empty group, that is match by position
	if (pEffectInfo->effectName() == "BasicRate") {
	//	LOGF(Priority::INFO, "'%s' in gmm group '%s'",
	//			pEffectInfo->effectName().c_str(), "");
		return "";
	}
	map<const string, const string>::const_iterator it =
		EffectFactory::GMM_GROUPS.find(pEffectInfo->effectName());
	// if group is not specified, map to base name
	if (it == EffectFactory::GMM_GROUPS.end()) {
	//	LOGF(Priority::INFO, "'%s' in gmm group '%s'",
	//			pEffectInfo->effectName().c_str(), pEffectInfo->effectName().c_str());
		return pEffectInfo->effectName();
	}
	// otherwise return the group
	//LOGF(Priority::INFO, "'%s' in gmm group '%s'",
	//		pEffectInfo->effectName().c_str(), it->second.c_str());
	return it->second;
}
const map<const string, const string> EffectFactory::GMM_GROUPS = init_groups();
map<const string, const string> EffectFactory::init_groups() {
	map<const string, const string> map;
	map.insert(make_pair("recip", "recip"));
	map.insert(make_pair("newrecip", "recip"));
	map.insert(make_pair("realrecip", "recip"));
	map.insert(make_pair("persistrecip", "recip"));
	//
	map.insert(make_pair("transTrip", "transTrip"));
	//
	map.insert(make_pair("egoX", "egoX"));
	map.insert(make_pair("egoX_gmm", "egoX"));
	map.insert(make_pair("outdeg", "egoX"));
	//
	map.insert(make_pair("simX", "simX"));
	map.insert(make_pair("simX_gmm", "simX"));
	map.insert(make_pair("totSim", "simX"));
	map.insert(make_pair("totSim_gmm", "simX"));
	map.insert(make_pair("avSim_gmm", "avSim"));
	map.insert(make_pair("avAlt_gmm", "avAlt"));
	map.insert(make_pair("totAlt_gmm", "totAlt"));
	map.insert(make_pair("maxAlt_gmm", "maxAlt"));
	map.insert(make_pair("minAlt_gmm", "minAlt"));
	return map;
}

/**
 * Constructor.
 * @param[in] pData the data this factory will create effects for
 */
EffectFactory::EffectFactory(const Data * pData)
{
	this->lpData = pData;
}

/**
 * Creates and returns a concrete effect of the Effect class hierarchy
 * corresponding to the given generic effect descriptor.
 */
Effect * EffectFactory::createEffect(const EffectInfo * pEffectInfo) const
{
	Effect * pEffect = 0;
	string effectName = pEffectInfo->effectName();

	// Defined so we can later on differentiate between effects for
	// continuous and discrete dependent behavior variables
	string variableName = pEffectInfo->variableName();
	ContinuousLongitudinalData * pContinuousData =
        dynamic_cast<ContinuousLongitudinalData *>(this->lpData->pContinuousData(variableName));

	// Handle the user-defined interaction effects first.

	if (pEffectInfo->pEffectInfo1())
	{
		// The info object of the first interacting effect is defined,
		// which means that we have a user-defined interaction effect.

		Effect *pEffect1 = this->createEffect(pEffectInfo->pEffectInfo1());
		Effect *pEffect2 = this->createEffect(pEffectInfo->pEffectInfo2());
		Effect *pEffect3 = 0;
		if (pEffectInfo->pEffectInfo3())
		{
			pEffect3 = this->createEffect(pEffectInfo->pEffectInfo3());
		}


		NetworkEffect * pNetworkEffect1 =
			dynamic_cast<NetworkEffect *>(pEffect1);
		BehaviorEffect * pBehaviorEffect1 =
			dynamic_cast<BehaviorEffect *>(pEffect1);

		if (pNetworkEffect1)
		{
			NetworkEffect * pNetworkEffect2 =
				dynamic_cast<NetworkEffect *>(pEffect2);
			NetworkEffect * pNetworkEffect3 = 0;
			if (pEffect3)
			{
				 pNetworkEffect3 = dynamic_cast<NetworkEffect *>(pEffect3);
			}
			pEffect = new NetworkInteractionEffect(pEffectInfo,
				pNetworkEffect1,
				pNetworkEffect2,
				pNetworkEffect3);
		}
		else
		{
			BehaviorEffect * pBehaviorEffect2 =
				dynamic_cast<BehaviorEffect *>(pEffect2);
			BehaviorEffect * pBehaviorEffect3 = 0;
			if (pEffect3)
			{
				pBehaviorEffect3 = dynamic_cast<BehaviorEffect *>(pEffect3);
			}

			pEffect = new BehaviorInteractionEffect(pEffectInfo,
				pBehaviorEffect1,
				pBehaviorEffect2,
				pBehaviorEffect3);
		}
	}
	else if (effectName == "density")
	{
		pEffect = new DensityEffect(pEffectInfo);
	}
	else if (effectName == "recip")
	{
		pEffect = new ReciprocityEffect(pEffectInfo);
	}
	else if (effectName == "newrecip")
	{
		pEffect = new ReciprocityGMMEffect(pEffectInfo,
				ReciprocityGMMEffect::NEW);
	}
	else if (effectName == "realrecip")
	{
		pEffect = new ReciprocityGMMEffect(pEffectInfo,
				ReciprocityGMMEffect::REAL);
	}
	else if (effectName == "persistrecip")
	{
		pEffect = new ReciprocityGMMEffect(pEffectInfo,
				ReciprocityGMMEffect::PERSISTENT);
	}
	else if (effectName == "transTrip1")
	{
		pEffect = new TransitiveTripletsEffect(pEffectInfo,true,false);
	}
	else if (effectName == "transTrip2")
	{
		pEffect = new TransitiveTripletsEffect(pEffectInfo,false,true);
	}
	else if (effectName == "transTrip")
	{
		pEffect = new TransitiveTripletsEffect(pEffectInfo,true,true);
	}
	else if (effectName == "agreetrans")
	{
		pEffect = new AgreementTransitivityEffect(pEffectInfo);
	}
	else if (effectName == "realtrans")
	{
		pEffect = new RealTransitivityEffect(pEffectInfo);
	}
	else if (effectName == "transTriads")
	{
		pEffect = new TransitiveTriadsEffect(pEffectInfo);
	}
	else if (effectName == "transMedTrip")
	{
		pEffect = new TransitiveMediatedTripletsEffect(pEffectInfo);
	}
	else if (effectName == "transRecTrip")
	{
		pEffect = new TransitiveReciprocatedTripletsEffect(pEffectInfo);
	}
	else if (effectName == "transRecTrip2")
	{
		pEffect = new TransitiveReciprocatedTriplets2Effect(pEffectInfo);
	}
	else if (effectName == "cycle3")
	{
		pEffect = new ThreeCyclesEffect(pEffectInfo);
	}
	else if (effectName == "transTies")
	{
		pEffect = new TransitiveTiesEffect(pEffectInfo);
	}
	else if (effectName == "between")
	{
		pEffect = new BetweennessEffect(pEffectInfo);
	}
	else if (effectName == "balance")
	{
		pEffect = new BalanceEffect(pEffectInfo);
	}
	else if (effectName == "nbrDist2")
	{
		pEffect = new DistanceTwoEffect(pEffectInfo, 1);
	}
	else if (effectName == "nbrDist2twice")
	{
		pEffect = new DistanceTwoEffect(pEffectInfo, 2);
	}
	else if (effectName == "denseTriads")
	{
		pEffect = new DenseTriadsEffect(pEffectInfo);
	}
	else if (effectName == "inPop")
	{
		pEffect = new IndegreePopularityEffect(pEffectInfo, false, false);
	}
	else if (effectName == "inPop.c")
	{
		pEffect = new IndegreePopularityEffect(pEffectInfo, false, true);
	}
	else if (effectName == "inPopSqrt")
	{
		pEffect = new IndegreePopularityEffect(pEffectInfo, true, false);
	}
	else if (effectName == "outPop")
	{
		pEffect = new OutdegreePopularityEffect(pEffectInfo, false, false);
	}
	else if (effectName == "outPop.c")
	{
		pEffect = new OutdegreePopularityEffect(pEffectInfo, false, true);
	}
	else if (effectName == "outPopSqrt")
	{
		pEffect = new OutdegreePopularityEffect(pEffectInfo, true, false);
	}
	else if (effectName == "reciPop")
	{
		pEffect = new RecipdegreePopularityEffect(pEffectInfo, false);
	}
	else if (effectName == "reciPopSqrt")
	{
		pEffect = new RecipdegreePopularityEffect(pEffectInfo, true);
	}
	else if (effectName == "inAct")
	{
		pEffect = new IndegreeActivityEffect(pEffectInfo, false, false);
	}
	else if (effectName == "inAct.c")
	{
		pEffect = new IndegreeActivityEffect(pEffectInfo, false, true);
	}
	else if (effectName == "inActSqrt")
	{
		pEffect = new IndegreeActivityEffect(pEffectInfo, true, false);
	}
	else if (effectName == "outAct")
	{
		pEffect = new OutdegreeActivityEffect(pEffectInfo, false);
	}
	else if (effectName == "outAct.c")
	{
		pEffect = new OutdegreeActivityEffect(pEffectInfo, true);
	}
	else if (effectName == "outActSqrt")
	{
		pEffect = new OutdegreeActivitySqrtEffect(pEffectInfo);
	}
	else if (effectName == "reciAct")
	{
		pEffect = new RecipdegreeActivityEffect(pEffectInfo);
	}
	else if (effectName == "degPlus")
	{
		pEffect = new BothDegreesEffect(pEffectInfo, false);
	}
	else if (effectName == "degPlus.c")
	{
		pEffect = new BothDegreesEffect(pEffectInfo, true);
	}
	else if (effectName == "avDeg")
	{
		pEffect = new AverageDegreeEffect(pEffectInfo);
	}
	else if (effectName == "outTrunc")
	{
		pEffect = new TruncatedOutdegreeEffect(pEffectInfo, true, false);
	}
	else if (effectName == "outTrunc2")
	{
		pEffect = new TruncatedOutdegreeEffect(pEffectInfo, true, false);
	}
	else if (effectName == "outMore")
	{
		pEffect = new TruncatedOutdegreeEffect(pEffectInfo, false, false);
	}
	else if (effectName == "outIso")
	{
		pEffect = new TruncatedOutdegreeEffect(pEffectInfo, true, true);
	}
	else if (effectName == "outInv")
	{
		pEffect = new InverseOutdegreeEffect(pEffectInfo);
	}
	else if (effectName == "outSqInv")
	{
		pEffect = new InverseSquaredOutdegreeEffect(pEffectInfo);
	}
	else if (effectName == "outOutAss")
	{
		pEffect = new OutOutDegreeAssortativityEffect(pEffectInfo);
	}
	else if (effectName == "outInAss")
	{
		pEffect = new OutInDegreeAssortativityEffect(pEffectInfo);
	}
	else if (effectName == "inOutAss")
	{
		pEffect = new InOutDegreeAssortativityEffect(pEffectInfo);
	}
	else if (effectName == "inInAss")
	{
		pEffect = new InInDegreeAssortativityEffect(pEffectInfo);
	}
	else if (effectName == "outAct_ego")
	{
		string networkName = pEffectInfo->variableName();
		AlterFunction * pFirstFunction =
			new EgoOutDegreeFunction(pEffectInfo->variableName());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->variableName(),
				AVERAGE_OUT_DEGREE);
		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "inAct_ego")
	{
		string networkName = pEffectInfo->variableName();
		AlterFunction * pFirstFunction =
			new EgoInDegreeFunction(pEffectInfo->variableName());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->variableName(),
				AVERAGE_IN_DEGREE);
		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "reciAct_ego")
	{
		string networkName = pEffectInfo->variableName();
		AlterFunction * pFirstFunction =
			new EgoRecipDegreeFunction(pEffectInfo->variableName());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->variableName(),
				AVERAGE_RECIP_DEGREE);
		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "X")
	{
		pEffect = new DyadicCovariateMainEffect(pEffectInfo);
	}
	else if (effectName == "XRecip")
	{
		pEffect = new DyadicCovariateReciprocityEffect(pEffectInfo);
	}
	else if (effectName == "WWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo, true, true);
	}
	else if (effectName == "cyWWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo, false, false);
	}
	else if (effectName == "InWWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo, false, true);
	}
	else if (effectName == "OutWWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo, true, false);
	}
	else if (effectName == "WXX")
	{
		pEffect = new WXXClosureEffect(pEffectInfo);
	}
	else if (effectName == "XWX")
	{
		pEffect = new XWXClosureEffect(pEffectInfo, true, true);
	}
	else if (effectName == "XWX1")
	{
		pEffect = new XWXClosureEffect(pEffectInfo, true, false);
	}
	else if (effectName == "XWX2")
	{
		pEffect = new XWXClosureEffect(pEffectInfo, false, true);
	}
	else if (effectName == "XXW")
	{
		pEffect = new XXWClosureEffect(pEffectInfo, true, true);
	}
	else if (effectName == "altX")
	{
		pEffect = new CovariateAlterEffect(pEffectInfo, false, false, false);
	}
	else if (effectName == "altSqX")
	{
		pEffect = new CovariateAlterEffect(pEffectInfo, false, false, true);
	}
	else if (effectName == "altLThresholdX")
	{
		pEffect = new CovariateAlterEffect(pEffectInfo, true, false, false);
	}
	else if (effectName == "altRThresholdX")
	{
		pEffect = new CovariateAlterEffect(pEffectInfo, false, true, false);
	}
	else if (effectName == "egoX")
	{
		pEffect = new CovariateEgoEffect(pEffectInfo, false, false);
	}
	else if (effectName == "egoSqX")
	{
		pEffect = new CovariateEgoSquaredEffect(pEffectInfo);
	}
	else if (effectName == "egoX_gmm")
	{
		pEffect = new CovariateEgoEffect(pEffectInfo, false, false, true);
	}
	else if (effectName == "egoPlusAltX")
	{
		pEffect = new CovariateDiffEffect(pEffectInfo, false, 0);
	}
	else if (effectName == "egoPlusAltSqX")
	{
		pEffect = new CovariateDiffEffect(pEffectInfo, false, 2);
	}
	else if (effectName == "egoLThresholdX")
	{
		pEffect = new CovariateEgoEffect(pEffectInfo, true, false);
	}
	else if (effectName == "egoRThresholdX")
	{
		pEffect = new CovariateEgoEffect(pEffectInfo, false, true);
	}
	else if (effectName == "degAbsDiffX")
	{
		pEffect = new CovariateEgoDiffEffect(pEffectInfo, true, true);
	}
	else if (effectName == "degPosDiffX")
	{
		pEffect = new CovariateEgoDiffEffect(pEffectInfo, true, false);
	}
	else if (effectName == "degNegDiffX")
	{
		pEffect = new CovariateEgoDiffEffect(pEffectInfo, false, true);
	}
	else if (effectName == "diffX")
	{
		pEffect = new CovariateDiffEffect(pEffectInfo, true, 0);
	}
	else if (effectName == "absDiffX")
	{
		pEffect = new CovariateDiffEffect(pEffectInfo, true, 1);
	}
	else if (effectName == "diffSqX")
	{
		pEffect = new CovariateDiffEffect(pEffectInfo, true, 2);
	}
	else if (effectName == "egoDiffX")
	{
		pEffect = new CovariateDiffEgoEffect(pEffectInfo);
	}
	else if (effectName == "simX")
	{
		pEffect = new CovariateSimilarityEffect(pEffectInfo, false);
	}
	else if (effectName == "simX_gmm")
	{
		pEffect = new CovariateSimilarityEffect(pEffectInfo, false, true);
	}
	else if (effectName == "simRecipX")
	{
		pEffect = new CovariateSimilarityEffect(pEffectInfo, true);
	}
	else if (effectName == "simRecipX_gmm")
	{
		pEffect = new CovariateSimilarityEffect(pEffectInfo, true, true);
	}
	else if (effectName == "simXTransTrip")
	{
		pEffect = new SimilarityTransitiveTripletsEffect(pEffectInfo, false);
	}
	else if (effectName == "sameX")
	{
		pEffect = new SameCovariateEffect(pEffectInfo, false);
	}
	else if (effectName == "higher")
	{
		pEffect = new HigherCovariateEffect(pEffectInfo);
	}
	else if (effectName == "sameXRecip")
	{
		pEffect = new SameCovariateEffect(pEffectInfo, true);
	}
	else if (effectName == "sameXTransTrip")
	{
		pEffect = new SameCovariateTransitiveTripletsEffect(pEffectInfo, true);
	}
	else if (effectName == "diffXTransTrip")
	{
		pEffect = new SameCovariateTransitiveTripletsEffect(pEffectInfo, false);
	}
	else if (effectName == "sameXTransRecTrip")
	{
		pEffect = new SameCovariateTransitiveReciprocatedTripletsEffect(pEffectInfo, true);
	}
	else if (effectName == "diffXTransRecTrip")
	{
		pEffect = new SameCovariateTransitiveReciprocatedTripletsEffect(pEffectInfo, false);
	}
	else if (effectName == "homXTransRecTrip")
	{
		pEffect = new HomCovariateTransitiveTripletsEffect(pEffectInfo, true);
	}
	else if (effectName == "inPopX")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new CovariateDegreeFunction(networkName, covariateName,
					false, true, false, (pEffectInfo->internalEffectParameter()>=2)),
			new CovariateDegreeFunction(networkName, covariateName,
					true, true, false, (pEffectInfo->internalEffectParameter()>=2)));
	}
	else if (effectName == "outPopX")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new CovariateDegreeFunction(networkName, covariateName,
					false, false, false, (pEffectInfo->internalEffectParameter()>=2)),
			new CovariateDegreeFunction(networkName, covariateName,
					true, false, false, (pEffectInfo->internalEffectParameter()>=2)));
	}
	else if (effectName == "inActX")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new CovariateDegreeFunction(networkName, covariateName,
					false, true, true, (pEffectInfo->internalEffectParameter()>=2)),
			new CovariateDegreeFunction(networkName, covariateName,
					true, true, true, (pEffectInfo->internalEffectParameter()>=2)));
	}
	else if (effectName == "outActX")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new CovariateDegreeFunction(networkName, covariateName,
					false, false, true, (pEffectInfo->internalEffectParameter()>=2)),
			new CovariateDegreeFunction(networkName, covariateName,
					true, false, true, (pEffectInfo->internalEffectParameter()>=2)));
	}
	else if (effectName == "sameXInPop")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new SameCovariateInTiesFunction(networkName, covariateName, true, true, false);
		AlterFunction * pStatisticFunction =
			new SameCovariateInTiesFunction(networkName, covariateName, true, true, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "diffXInPop")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new SameCovariateInTiesFunction(networkName, covariateName, false, true, false);
		AlterFunction * pStatisticFunction =
			new SameCovariateInTiesFunction(networkName, covariateName, false, true, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "sameXOutAct")
	{
		pEffect = new SameCovariateActivityEffect(pEffectInfo, true, false);
	}
	else if (effectName == "diffXOutAct")
	{
		pEffect = new SameCovariateActivityEffect(pEffectInfo, false, false);
	}
	else if (effectName == "homXOutAct")
	{
		pEffect = new HomCovariateActivityEffect(pEffectInfo, true);
	}
	else if (effectName == "altXOutAct")
	{
		pEffect = new AlterCovariateActivityEffect(pEffectInfo);
	}
	else if (effectName == "sameXReciAct")
	{
		pEffect = new SameCovariateActivityEffect(pEffectInfo, true, true);
	}
	else if (effectName == "diffXReciAct")
	{
		pEffect = new SameCovariateActivityEffect(pEffectInfo, false, true);
	}
	else if (effectName == "transTripX")
	{
		pEffect = new CovariateTransitiveTripletsEffect(pEffectInfo);
	}
	else if (effectName == "homXTransTrip")
	{
		pEffect = new HomCovariateTransitiveTripletsEffect(pEffectInfo, false);
	}
	else if (effectName == "jumpXTransTrip")
	{
		pEffect = new JumpCovariateTransitiveTripletsEffect(pEffectInfo, false);
	}
	else if (effectName == "egoXaltX")
	{
		pEffect = new CovariateEgoAlterEffect(pEffectInfo, false);
	}
	else if (effectName == "egoXaltXRecip")
	{
		pEffect = new CovariateEgoAlterEffect(pEffectInfo, true);
	}
	else if (effectName == "IndTies")
	{
		pEffect = new CovariateIndirectTiesEffect(pEffectInfo);
	}
	else if (effectName == "sameXCycle4")
	{
		pEffect = new SameCovariateFourCyclesEffect(pEffectInfo, true);
	}
	else if (effectName == "avGroupEgoX")
	{
		pEffect = new AverageGroupEgoEffect(pEffectInfo);
	}
	else if (effectName == "cycle4")
	{
		pEffect = new FourCyclesEffect(pEffectInfo, true);
	}
	else if (effectName == "sharedPop")
	{
		pEffect = new FourCyclesEffect(pEffectInfo, false);
	}
	else if (effectName == "cycle4ND")
	{
		pEffect = new FourCyclesEffect(pEffectInfo, false);
	}
	else if ((effectName == "gwespFF")||(effectName == "gwesp"))
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pTwoPathTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespFB")
 	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pInStarTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());

		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespBF")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pOutStarTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespBB")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pReverseTwoPathTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespRR")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pRRTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwdspFF")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pTwoPathTable;
		pEffect = new GwdspEffect(pEffectInfo, mytable, true);
	}
	else if (effectName == "gwdspFB")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pInStarTable;
		pEffect = new GwdspEffect(pEffectInfo, mytable, false);
	}
	else if (effectName == "inStructEq")
	{
		pEffect = new InStructuralEquivalenceEffect(pEffectInfo);
	}
	else if (effectName == "Jin")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new InJaccardFunction(pEffectInfo->variableName()));
	}
	else if (effectName == "Jout")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new OutJaccardFunction(pEffectInfo->variableName()));
	}
	else if (effectName == "crprod")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new OutTieFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "crprodRecip")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new InTieFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "crprodMutual")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ProductFunction(
				new OutTieFunction(pEffectInfo->interactionName1()),
				new InTieFunction(pEffectInfo->interactionName1())));
	}
	else if (effectName == "inPopIntn")
	{
		AlterFunction * pFirstFunction =
			new InDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_IN_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "inActIntn")
	{
		AlterFunction * pFirstFunction =
			new EgoInDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_IN_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
    else if (effectName == "crprodInActIntn")
    {
        AlterFunction * pFirstFunction =
            new EgoInDegreeFunction(pEffectInfo->interactionName1());
        ConstantFunction * pSecondFunction =
            new ConstantFunction(pEffectInfo->interactionName1(),
                AVERAGE_IN_DEGREE);

        if (pEffectInfo->internalEffectParameter() == 2)
        {
            pFirstFunction = new IntSqrtFunction(pFirstFunction);
            pSecondFunction->pFunction(sqrt);
        }

        pEffect = new GenericNetworkEffect(pEffectInfo,
            new ProductFunction(
                new OutTieFunction(pEffectInfo->interactionName1()),
                new DifferenceFunction(pFirstFunction, pSecondFunction)));
    }
	else if (effectName == "outPopIntn")
	{
		AlterFunction * pFirstFunction =
			new OutDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_OUT_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "outActIntn")
	{
		AlterFunction * pFirstFunction =
			new EgoOutDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_OUT_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "absOutDiffIntn")
	{
		AlterFunction * pAlterOutdegreeFunction =
			new OutDegreeFunction(pEffectInfo->interactionName1());
		AlterFunction * pEgoOutdegreeFunction =
			new EgoOutDegreeFunction(pEffectInfo->interactionName1());
		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pAlterOutdegreeFunction =
				new IntSqrtFunction(pAlterOutdegreeFunction);
			pEgoOutdegreeFunction =
				new IntSqrtFunction(pEgoOutdegreeFunction);
		}
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new AbsDiffFunction(pAlterOutdegreeFunction,pEgoOutdegreeFunction));
	}
	else if (effectName == "avDegIntn")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
						new DegreeFunction(pEffectInfo->interactionName1(), 
									pEffectInfo->internalEffectParameter()));
	}
	else if (effectName == "doubleOutAct")
	{
		string firstNetworkName = pEffectInfo->variableName();
		string secondNetworkName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new DoubleOutActFunction(firstNetworkName, secondNetworkName,
					pEffectInfo->internalEffectParameter(), true);
		AlterFunction * pStatisticFunction =
			new DoubleOutActFunction(firstNetworkName, secondNetworkName,
					pEffectInfo->internalEffectParameter(), false);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
// doubleOutAct and doubleInPop are implemented differently
// just to explore the differences.
	else if (effectName == "doubleInPop")
	{
		pEffect = new DoubleInPopEffect(pEffectInfo);
	}
	else if (effectName == "JinMix")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new InJaccardFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "JoutMix")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new OutJaccardFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "gwespFFMix")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pTwoPathTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo-> interactionName1(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if ((effectName == "gwespFBMix")||(effectName == "gwespMix"))
 	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pInStarTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->interactionName1(),
				mytable, pEffectInfo->internalEffectParameter());
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespBFMix")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pOutStarTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->interactionName1(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespBBMix")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pReverseTwoPathTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->interactionName1(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespRRMix")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pRRTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->interactionName1(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "both")
	{
		AlterFunction * pAlterIndegreeFunction =
			new InDegreeFunction(pEffectInfo->interactionName1());
		AlterFunction * pEgoIndegreeFunction =
			new EgoInDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pFirstConstantFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_IN_DEGREE);
		ConstantFunction * pSecondConstantFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_IN_DEGREE);
		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pAlterIndegreeFunction =
				new IntSqrtFunction(pAlterIndegreeFunction);
			pEgoIndegreeFunction =
				new IntSqrtFunction(pEgoIndegreeFunction);
			pFirstConstantFunction->pFunction(sqrt);
			pSecondConstantFunction->pFunction(sqrt);
		}
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ProductFunction(
				new DifferenceFunction(pAlterIndegreeFunction,
					pFirstConstantFunction),
				new DifferenceFunction(pEgoIndegreeFunction,
					pSecondConstantFunction)));
	}
	else if (effectName == "betweenPop")
	{
		AlterFunction * pFunction =
			new BetweennessFunction(pEffectInfo->interactionName1());

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFunction = new IntSqrtFunction(pFunction);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo, pFunction);
	}
	else if (effectName == "from")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new InStarFunction(pEffectInfo->interactionName1(),
							(pEffectInfo->internalEffectParameter() >= 2)));
	}
	else if (effectName == "fromMutual")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ReciprocatedTwoPathFunction(pEffectInfo->interactionName1(),
							(pEffectInfo->internalEffectParameter() >= 2)));
	}
	else if (effectName == "to")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
				new MixedTwoStepFunction(
						pEffectInfo->interactionName1(),
						pEffectInfo->variableName(),
						FORWARD, FORWARD, pEffectInfo->internalEffectParameter()));
	}
	else if (effectName == "toBack") // formerly mixedInWX
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
				new MixedTwoStepFunction(
						pEffectInfo->interactionName1(),
						pEffectInfo->variableName(),
						BACKWARD, FORWARD, pEffectInfo->internalEffectParameter()));
		}
	else if (effectName == "toRecip")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
				new MixedTwoStepFunction(
						pEffectInfo->interactionName1(),
						pEffectInfo->variableName(),
						RECIPROCAL, FORWARD, pEffectInfo->internalEffectParameter()));
	}
	else if (effectName == "toAny")
	{
		pEffect = new MixedOnlyTwoPathEffect(pEffectInfo);
	}
	else if (effectName == "cl.XWX")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new SumFunction(
			new MixedTwoStepFunction(
						pEffectInfo->variableName(),
						pEffectInfo->interactionName1(),
						FORWARD, FORWARD, 0),
			new MixedTwoStepFunction(
						pEffectInfo->variableName(),
						pEffectInfo->interactionName1(),
						FORWARD, BACKWARD, 0)));
	}
	else if (effectName == "cl.XWX1")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
					new MixedTwoStepFunction(
						pEffectInfo->variableName(),
						pEffectInfo->interactionName1(),
						FORWARD, FORWARD, 0));
	}
	else if (effectName == "cl.XWX2")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new  MixedTwoStepFunction(
						pEffectInfo->variableName(),
						pEffectInfo->interactionName1(),
						FORWARD, BACKWARD, 0));
	}
	else if (effectName == "mixedInXW")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
				new MixedTwoStepFunction(
						pEffectInfo->variableName(),
						pEffectInfo->interactionName1(),
						BACKWARD, FORWARD, pEffectInfo->internalEffectParameter()));
	}
	/*
	 * Mixed two step effects with flexible iterators (added by Christoph)
	 */
	else if (effectName == "transTrip.FR")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new TwoStepFunction(pEffectInfo->variableName(),
					FORWARD, RECIPROCAL));
	}
	else if (effectName == "transTrip.EE")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoStepFunction(pEffectInfo->variableName(), // not mixed, TwoStepFunction equivalent not implemented
					pEffectInfo->variableName(),
					EITHER, EITHER, 0));
	}
	else if (effectName == "transTrip.FE")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoStepFunction(pEffectInfo->variableName(), // not mixed, TwoStepFunction equivalent not implemented
					pEffectInfo->variableName(),
					FORWARD, EITHER, 0));
	}
	else if (effectName == "WWX.EE")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoStepFunction(pEffectInfo->interactionName1(), // not mixed
					pEffectInfo->interactionName1(),
					EITHER, EITHER, 0));
	}
	else if (effectName == "WXX.FE")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoStepFunction(pEffectInfo->interactionName1(),
					pEffectInfo->variableName(),
					FORWARD, EITHER, 0));
	}
	else if (effectName == "XWX.ER")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoStepFunction(pEffectInfo->variableName(),
					pEffectInfo->interactionName1(),
					EITHER, RECIPROCAL, 0));
	}
	else if (effectName == "WWX.FR")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoStepFunction(pEffectInfo->interactionName1(),  // not mixed
					pEffectInfo->interactionName1(),
					FORWARD, RECIPROCAL, 0));
	}
	else if (effectName == "WXX.ER")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoStepFunction(pEffectInfo->interactionName1(),
					pEffectInfo->variableName(),
					EITHER, RECIPROCAL, 0));
	}
	else if (effectName == "XWX.FE")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoStepFunction(pEffectInfo->variableName(),
					pEffectInfo->interactionName1(),
					FORWARD, EITHER, 0));
	}
	else if (effectName == "outOutActIntn")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new OutActDistance2Function(pEffectInfo->interactionName1(),
							pEffectInfo->variableName(),
							pEffectInfo->internalEffectParameter(), false, false, false));
	}
	else if (effectName == "outOutAvIntn")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new OutActDistance2Function(pEffectInfo->interactionName1(),
							pEffectInfo->variableName(),
							pEffectInfo->internalEffectParameter(), false, false, true));
	}
	else if (effectName == "sharedTo")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
						new MixedThreeCyclesFunction(pEffectInfo->interactionName1(),
							pEffectInfo->variableName(),
							pEffectInfo->internalEffectParameter()));
	}
	else if (effectName == "toU")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new WeightedMixedTwoPathFunction(pEffectInfo->interactionName1(),
					pEffectInfo->variableName(),
					pEffectInfo->interactionName2(), false),
			new WeightedMixedTwoPathFunction(pEffectInfo->interactionName1(),
					pEffectInfo->variableName(),
					pEffectInfo->interactionName2(), true));
	}
	else if (effectName == "covNetNet")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new InStarFunction(networkName, false),
				0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new InStarFunction(networkName, false),
					0)));
	}
	else if (effectName == "homCovNetNet")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new SameCovariateInStarFunction(networkName,
												covariateName, false), 0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new SameCovariateInStarFunction(networkName,
												covariateName, true), 0)));
	}
	else if (effectName == "jumpWWClosure")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				0,
				new SameCovariateTwoPathFunction(networkName,
										covariateName, false)),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					0,
					new SameCovariateTwoPathFunction(networkName,
										covariateName, true))));
	}
	else if (effectName == "jumpFrom")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				0,
				new SameCovariateInStarFunction(networkName,
										covariateName, false)),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					0,
					new SameCovariateInStarFunction(networkName,
										covariateName, true))));
	}
	else if (effectName == "contrastCovNetNet")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new DifferentCovariateInStarFunction(networkName,
										covariateName, false, false),
				0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new DifferentCovariateInStarFunction(networkName,
										covariateName, true, false),
					0)));
	}
	else if (effectName == "allDifCovNetNet")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new DifferentCovariateInStarFunction(networkName,
										covariateName, false, true),
				0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new DifferentCovariateInStarFunction(networkName,
										covariateName, true, true),
					0)));
	}
	else if (effectName == "jumpSharedIn")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				0,
				new SameCovariateOutStarFunction(networkName,
										covariateName, false)),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					0,
					new SameCovariateOutStarFunction(networkName,
										covariateName, true))));
	}
	else if (effectName == "covNetNetIn")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new OutStarFunction(networkName),
				0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new OutStarFunction(networkName),
					0)));
	}
	else if (effectName == "homCovNetNetIn")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new SameCovariateOutStarFunction(networkName,
												covariateName, false), 0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new SameCovariateOutStarFunction(networkName,
												covariateName, true), 0)));
	}
	else if (effectName == "contrastCovNetNetIn")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new DifferentCovariateOutStarFunction(networkName,
										covariateName, false, false),
				0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new DifferentCovariateOutStarFunction(networkName,
										covariateName, true, false),
					0)));
	}
	else if (effectName == "allDifCovNetNetIn")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new DifferentCovariateOutStarFunction(networkName,
										covariateName, false, false),
				0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new DifferentCovariateOutStarFunction(networkName,
										covariateName, true, false),
					0)));
	}
	else if (effectName == "sameWXClosure")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new SameCovariateMixedTwoPathFunction(
							pEffectInfo->variableName(),
							networkName, covariateName, false),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new SameCovariateMixedTwoPathFunction(
							pEffectInfo->variableName(),
							networkName, covariateName, true)));
	}
	else if (effectName == "homWXClosure")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new SameCovariateMixedTwoPathFunction(
								pEffectInfo->variableName(),
								networkName, covariateName, false), 0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new SameCovariateMixedTwoPathFunction(
							pEffectInfo->variableName(),
							networkName, covariateName, true), 0)));
	}
	else if (effectName == "jumpWXClosure")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				0,
				new SameCovariateMixedTwoPathFunction(
								pEffectInfo->variableName(),
								networkName, covariateName, false)),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					0,
					new SameCovariateMixedTwoPathFunction(
							pEffectInfo->variableName(),
							networkName, covariateName, true))));
	}
	else if (effectName == "closure")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new TwoPathFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "cyClosure")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ReverseTwoPathFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "sharedIn")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new OutStarFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "from.w.ind")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new InStarsTimesDegreesFunction(pEffectInfo->interactionName1(),
						pEffectInfo->interactionName2(),
						pEffectInfo->internalEffectParameter()));
	}
	else if (effectName == "linear")
	{
		pEffect = new LinearShapeEffect(pEffectInfo);
	}
	else if (effectName == "quad")
	{
		pEffect = new QuadraticShapeEffect(pEffectInfo);
	}
	else if (effectName == "constant")
	{
		pEffect = new ConstantEffect(pEffectInfo);
	}
	else if ((effectName == "threshold") || (effectName == "threshold2") ||
	         (effectName == "threshold3") || (effectName == "threshold4"))
	{
		pEffect = new ThresholdShapeEffect(pEffectInfo);
	}
	else if (effectName == "simAllNear")
	{
		pEffect = new AllSimilarityEffect(pEffectInfo, true);
	}
	else if (effectName == "simAllFar")
	{
		pEffect = new AllSimilarityEffect(pEffectInfo, false);
	}
	else if (effectName == "avSim")
	{
		pEffect = new SimilarityEffect(pEffectInfo, true, false, false, true, true);
	}
	else if (effectName == "avAttHigher")
	{
		pEffect = new SimilarityEffect(pEffectInfo, true, false, false, true, false);
	}
	else if (effectName == "avAttLower")
	{
		pEffect = new SimilarityEffect(pEffectInfo, true, false, false, false, true);
	}
	else if (effectName == "avSim_gmm")
	{
		pEffect = new SimilarityEffect(pEffectInfo, true, false, false, true, true, true);
	}
	else if (effectName == "totSim")
	{
		pEffect = new SimilarityEffect(pEffectInfo, false, false, false, true, true);
	}
	else if (effectName == "totAttHigher")
	{
		pEffect = new SimilarityEffect(pEffectInfo, false, false, false, true, false);
	}
	else if (effectName == "totAttLower")
	{
		pEffect = new SimilarityEffect(pEffectInfo, false, false, false, false, true);
	}
	else if (effectName == "avInSim")
	{
		pEffect = new SimilarityIndegreeEffect(pEffectInfo, true, false);
	}
	else if (effectName == "totInSim")
	{
		pEffect = new SimilarityIndegreeEffect(pEffectInfo, false, false);
	}
	else if (effectName == "indeg")
	{
		if (pContinuousData)
		{
			pEffect = new IndegreeContinuousEffect(pEffectInfo, false);
		}
		else
		{
		pEffect = new IndegreeEffect(pEffectInfo);
		}
	}
	else if (effectName == "indegSqrt")
	{
		if (pContinuousData)
		{
			pEffect = new IndegreeContinuousEffect(pEffectInfo, true);
		}
	}
	else if (effectName == "outdeg")
	{
		if (pContinuousData)
		{
			pEffect = new OutdegreeContinuousEffect(pEffectInfo, false);
		}
		else
		{
		pEffect = new OutdegreeEffect(pEffectInfo);
		}
	}
	else if (effectName == "outdegSqrt")
	{
		if (pContinuousData)
		{
			pEffect = new OutdegreeContinuousEffect(pEffectInfo, true);
		}
	}
	else if (effectName == "isolate")
	{
		pEffect = new IsolateEffect(pEffectInfo, true);
	}
	else if (effectName == "outIsolate")
	{
		pEffect = new IsolateEffect(pEffectInfo, false);
	}
	else if (effectName == "isolateNet")
	{
		pEffect = new IsolateNetEffect(pEffectInfo);
	}
	else if (effectName == "antiIso")
	{
		pEffect = new AntiIsolateEffect(pEffectInfo, true, 1);
	}
	else if (effectName == "antiInIso")
	{
		pEffect = new AntiIsolateEffect(pEffectInfo, false, 1);
	}
	else if ((effectName == "antiInIso2") || (effectName == "in2Plus"))
	{
		pEffect = new AntiIsolateEffect(pEffectInfo, false, 2);
	}
	else if (effectName == "in3Plus")
	{
		pEffect = new AntiIsolateEffect(pEffectInfo, false, 3);
	}
	else if (effectName == "isolateOut")
	{
		pEffect = new IsolateOutContinuousEffect(pEffectInfo);
	}
	else if (effectName == "isolatePop")
	{
		pEffect = new IsolatePopEffect(pEffectInfo, true);
	}
	else if (effectName == "inIsDegree")
	{
		pEffect = new InIsolateDegreeEffect(pEffectInfo);
	}
	else if (effectName == "nonPrimary")
	{
		pEffect = new PrimaryCompressionEffect(pEffectInfo, false, false);
	}
	else if (effectName == "nonPCompress")
	{
		pEffect = new PrimaryCompressionEffect(pEffectInfo, false, true);
	}
	else if (effectName == "primary")
	{
		pEffect = new PrimaryCompressionEffect(pEffectInfo, true, false);
	}
	else if (effectName == "primCompress")
	{
		pEffect = new PrimaryCompressionEffect(pEffectInfo, true, true);
	}
	else if (effectName == "primDegAct")
	{
		pEffect = new PrimarySettingEffect(pEffectInfo, false, false, false, false);
	}
	else if (effectName == "primDegActDiff")
	{
		pEffect = new PrimarySettingEffect(pEffectInfo, true, false, false, false);
	}
	else if (effectName == "primDegActDiffSqrt")
	{
		pEffect = new PrimarySettingEffect(pEffectInfo, true, false, true, false);
	}
	else if (effectName == "primDegActSqrt")
	{
		pEffect = new PrimarySettingEffect(pEffectInfo, false, false, true, false);
	}
	else if (effectName == "primDegActLog")
	{
		pEffect = new PrimarySettingEffect(pEffectInfo, false, true, false, false);
	}
	else if (effectName == "primDegActInv")
	{
		pEffect = new PrimarySettingEffect(pEffectInfo, false, false, false, true);
	}
	else if (effectName == "settingSizeAct")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, false, false, false, false, false, false);
	}
	else if (effectName == "settingSizeActSqrt")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, false, false, true, false, false, false);
	}
	else if (effectName == "settingSizeActLog")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, false, true, false, false, false, false);
	}
	else if (effectName == "settingSizeActInv")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, false, false, false, true, false, false);
	}
	else if (effectName == "settingOppAct")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, true, false, false, false, false, false);
	}
	else if (effectName == "settingOppActSqrt")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, true, false, true, false, false, false);
	}
	else if (effectName == "settingOppActLog")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, true, true, false, false, false, false);
	}
	else if (effectName == "settingLogCreationAct")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, true, true, false, false, true, false);
	}
	else if (effectName == "settingOppActD")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, true, false, false, false, false, true);
	}
	else if (effectName == "settingOppActSqrtD")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, true, false, true, false, false, true);
	}
	else if (effectName == "settingOppActLogD")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, true, true, false, false, false, true);
	}
	else if (effectName == "settingLogCreationActD")
	{
		pEffect = new SettingSizeEffect(pEffectInfo, true, true, false, false, true, true);
	}
	else if (effectName == "avSimRecip")
	{
		pEffect = new ReciprocatedSimilarityEffect(pEffectInfo, true, false);
	}
	else if (effectName == "totSimRecip")
	{
		pEffect = new ReciprocatedSimilarityEffect(pEffectInfo, false, false);
	}
	else if (effectName == "avSimPopAlt")
	{
		pEffect = new SimilarityEffect(pEffectInfo, true, true, false, true, true);
	}
	else if (effectName == "totSimPopAlt")
	{
		pEffect = new SimilarityEffect(pEffectInfo, false, true, false, true, true);
	}
	else if (effectName == "avInSimPopAlt")
	{
		pEffect = new SimilarityIndegreeEffect(pEffectInfo, true, true);
	}
	else if (effectName == "totInSimPopAlt")
	{
		pEffect = new SimilarityIndegreeEffect(pEffectInfo, false, true);
	}
	else if (effectName == "popAlt")
	{
		pEffect = new PopularityAlterEffect(pEffectInfo);
	}
	else if (effectName == "avSimRecPop")
	{
		pEffect = new ReciprocatedSimilarityEffect(pEffectInfo, true, true);
	}
	else if (effectName == "totSimRecPop")
	{
		pEffect = new ReciprocatedSimilarityEffect(pEffectInfo, false, true);
	}
	else if (effectName == "avAlt")
	{
		if (pContinuousData)
		{
			pEffect = new AverageAlterContinuousEffect(pEffectInfo);
		}
		else
		{
		pEffect = new AverageAlterEffect(pEffectInfo, true, false);
	}
	}
	else if (effectName == "avAlt_gmm")
	{
		pEffect = new AverageAlterEffect(pEffectInfo, true, false, true);
	}
	else if (effectName == "avGroup")
	{
		pEffect = new AverageGroupEffect(pEffectInfo);
	}
	else if (effectName == "totAlt")
	{
		pEffect = new AverageAlterEffect(pEffectInfo, false, false);
	}
	else if (effectName == "totAlt_gmm")
	{
		pEffect = new AverageAlterEffect(pEffectInfo, false, false, true);
	}
	else if (effectName == "avAltPop")
	{
		pEffect = new AverageAlterEffect(pEffectInfo, true, true);
	}
	else if (effectName == "totAltPop")
	{
		pEffect = new AverageAlterEffect(pEffectInfo, false, true);
	}
	else if (effectName == "avRecAlt")
	{
		pEffect = new AverageReciprocatedAlterEffect(pEffectInfo, true);
	}
	else if (effectName == "totRecAlt")
	{
		pEffect = new AverageReciprocatedAlterEffect(pEffectInfo, false);
	}
	else if (effectName == "maxAlt")
	{
		if (pContinuousData)
		{
			pEffect = new MaxAlterContinuousEffect(pEffectInfo, false);
		}
		else
		{
		pEffect = new MaxAlterEffect(pEffectInfo, false);
	}
	}
	else if (effectName == "maxAlt_gmm")
	{
		pEffect = new MaxAlterEffect(pEffectInfo, false, true);
	}
	else if (effectName == "minAlt")
	{
		if (pContinuousData)
		{
			pEffect = new MaxAlterContinuousEffect(pEffectInfo, true);
		}
		else
		{
		pEffect = new MaxAlterEffect(pEffectInfo, true);
	}
	}
	else if (effectName == "minAlt_gmm")
	{
		pEffect = new MaxAlterEffect(pEffectInfo, true, true);
	}
	else if (effectName == "avInAlt")
	{
		pEffect = new AverageInAlterEffect(pEffectInfo, true);
	}
	else if (effectName == "totInAlt")
	{
		pEffect = new AverageInAlterEffect(pEffectInfo, false);
	}
	else if (effectName == "avAltDist2")
	{
		pEffect = new AverageAlterDist2Effect(pEffectInfo, true, true);
	}
	else if (effectName == "totAltDist2")
	{
		pEffect = new AverageAlterDist2Effect(pEffectInfo, false, false);
	}
	else if (effectName == "avTAltDist2")
	{
		pEffect = new AverageAlterDist2Effect(pEffectInfo, true, false);
	}
	else if (effectName == "totAAltDist2")
	{
		pEffect = new AverageAlterDist2Effect(pEffectInfo, false, true);
	}
	else if (effectName == "avInAltDist2")
	{
		pEffect = new AverageAlterInDist2Effect(pEffectInfo, true, true);
	}
	else if (effectName == "totInAltDist2")
	{
		pEffect = new AverageAlterInDist2Effect(pEffectInfo, false, false);
	}
	else if (effectName == "avTInAltDist2")
	{
		pEffect = new AverageAlterInDist2Effect(pEffectInfo, true, false);
	}
	else if (effectName == "totAInAltDist2")
	{
		pEffect = new AverageAlterInDist2Effect(pEffectInfo, false, true);
	}
	else if (effectName == "behDenseTriads")
	{
		pEffect = new DenseTriadsBehaviorEffect(pEffectInfo);
	}
	else if (effectName == "simDenseTriads")
	{
		pEffect = new DenseTriadsSimilarityEffect(pEffectInfo);
	}
	else if (effectName == "recipDeg")
	{
		if (pContinuousData)
		{
			pEffect = new ReciprocalDegreeContinuousEffect(pEffectInfo, true);
		}
		else
		{
		pEffect = new ReciprocalDegreeBehaviorEffect(pEffectInfo);
		}
	}
	else if (effectName == "nonrecipDeg")
	{
		pEffect = new ReciprocalDegreeContinuousEffect(pEffectInfo, false);
	}
	else if (effectName == "FFDeg")
	{
		pEffect = new DoubleDegreeBehaviorEffect(pEffectInfo, true, 0);
	}
	else if (effectName == "BBDeg")
	{
		pEffect = new DoubleDegreeBehaviorEffect(pEffectInfo, false, 1);
	}
	else if (effectName == "FBDeg")
	{
		pEffect = new DoubleDegreeBehaviorEffect(pEffectInfo, true, 1);
	}
	else if (effectName == "BFDeg")
	{
		pEffect = new DoubleDegreeBehaviorEffect(pEffectInfo, false, 0);
	}
	else if (effectName == "FRDeg")
	{
		pEffect = new DoubleDegreeBehaviorEffect(pEffectInfo, true, 2);
	}
	else if (effectName == "BRDeg")
	{
		pEffect = new DoubleDegreeBehaviorEffect(pEffectInfo, false, 2);
	}
//	else if (effectName == "RRDeg")
//	{
//		pEffect = new DoubleRecDegreeBehaviorEffect(pEffectInfo, 2); // leads to error
//	}
	else if (effectName == "avSimPopEgo")
	{
		pEffect = new SimilarityEffect(pEffectInfo, true, false, true, true, true);
	}
	else if (effectName == "effFrom")
	{
		if (pContinuousData)
		{
			pEffect = new MainCovariateContinuousEffect(pEffectInfo);
		}
		else
		{
		pEffect = new MainCovariateEffect(pEffectInfo);
		}
	}
	else if (effectName == "avSimEgoX")
	{
		pEffect = new InteractionCovariateEffect(pEffectInfo, true, false, false, false);
	}
	else if (effectName == "totSimEgoX")
	{
		pEffect = new InteractionCovariateEffect(pEffectInfo, false, true, false, false);
	}
	else if (effectName == "avAltEgoX")
	{
		pEffect = new InteractionCovariateEffect(pEffectInfo, false, false, true, false);
	}
	else if (effectName == "totAltEgoX")
	{
		pEffect = new InteractionCovariateEffect(pEffectInfo, false, false, false, true);
	}
	else if (effectName == "totSimAltX")
	{
		pEffect = new AltersCovariateTotSimEffect(pEffectInfo);
	}
	else if (effectName == "avSimAltX")
	{
		pEffect = new AltersCovariateAvSimEffect(pEffectInfo);
	}
	else if (effectName == "totAltAltX")
	{
		pEffect = new AltersCovariateAvAltEffect(pEffectInfo, false);
	}
	else if (effectName == "avAltAltX")
	{
		pEffect = new AltersCovariateAvAltEffect(pEffectInfo, true);
	}
	else if (effectName == "AltsAvAlt")
	{
		throw domain_error("Effect AltsAvAlt renamed to avXAlt");
	}
	else if (effectName == "minXAlt")
	{
		pEffect = new AltersCovariateMinimumEffect(pEffectInfo);
	}
	else if (effectName == "maxXAlt")
	{
		pEffect = new AltersCovariateMaximumEffect(pEffectInfo);
	}
	else if (effectName == "avXAlt")
	{
		pEffect = new AltersCovariateAverageEffect(pEffectInfo,true);
	}
	else if (effectName == "totXAlt")
	{
		pEffect = new AltersCovariateAverageEffect(pEffectInfo,false);
	}
	else if (effectName == "avXInAlt")
	{
		pEffect = new InAltersCovariateAverageEffect(pEffectInfo,true);
	}
	else if (effectName == "totXInAlt")
	{
		pEffect = new InAltersCovariateAverageEffect(pEffectInfo,false);
	}
	else if (effectName == "avXAltDist2")
	{
		pEffect = new AltersDist2CovariateAverageEffect(pEffectInfo,true,true);
	}
	else if (effectName == "totXAltDist2")
	{
		pEffect = new AltersDist2CovariateAverageEffect(pEffectInfo,false,false);
	}
	else if (effectName == "avTXAltDist2")
	{
		pEffect = new AltersDist2CovariateAverageEffect(pEffectInfo,true,false);
	}
	else if (effectName == "totAXAltDist2")
	{
		pEffect = new AltersDist2CovariateAverageEffect(pEffectInfo,false,true);
	}
	else if (effectName == "avXInAltDist2")
	{
		pEffect = new AltersInDist2CovariateAverageEffect(pEffectInfo,true,true);
	}
	else if (effectName == "totXInAltDist2")
	{
		pEffect = new AltersInDist2CovariateAverageEffect(pEffectInfo,false,false);
	}
	else if (effectName == "avTXInAltDist2")
	{
		pEffect = new AltersInDist2CovariateAverageEffect(pEffectInfo,true,false);
	}
	else if (effectName == "totAXInAltDist2")
	{
		pEffect = new AltersInDist2CovariateAverageEffect(pEffectInfo,false,true);
	}
	else if (effectName == "degAbsContrX")
	{
		pEffect = new CovariateContrastEffect(pEffectInfo, true, true);
	}
	else if (effectName == "degPosContrX")
	{
		pEffect = new CovariateContrastEffect(pEffectInfo, true, false);
	}
	else if (effectName == "degNegContrX")
	{
		pEffect = new CovariateContrastEffect(pEffectInfo, false, true);
	}
	else if (effectName == "avWAlt")
	{
		pEffect = new DyadicCovariateAvAltEffect(pEffectInfo,true, false, true);
	}
	else if (effectName == "avWInAlt")
	{
		pEffect = new DyadicCovariateAvAltEffect(pEffectInfo,true, false, false);
	}
	else if (effectName == "totWAlt")
	{
		pEffect = new DyadicCovariateAvAltEffect(pEffectInfo,false, false, true);
	}
	else if (effectName == "totWInAlt")
	{
		pEffect = new DyadicCovariateAvAltEffect(pEffectInfo,false, false, false);
	}
	else if (effectName == "avAltW")
	{
		pEffect = new DyadicCovariateAvAltEffect(pEffectInfo,true, true, true);
	}
	else if (effectName == "avInAltW")
	{
		pEffect = new DyadicCovariateAvAltEffect(pEffectInfo,true, true, false);
	}
	else if (effectName == "totAltW")
	{
		pEffect = new DyadicCovariateAvAltEffect(pEffectInfo,false, true, true);
	}
	else if (effectName == "totInAltW")
	{
		pEffect = new DyadicCovariateAvAltEffect(pEffectInfo,false, true, false);
	}
	else if (effectName == "avSimW")
	{
		pEffect = new SimilarityWEffect(pEffectInfo, true, false, false);
	}
	else if (effectName == "totSimW")
	{
		pEffect = new SimilarityWEffect(pEffectInfo, false, false, false);
	}
	else if (effectName == "altDist2")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		double parameter = pEffectInfo->internalEffectParameter();
		AlterFunction * pChangeFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, false, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, true, false);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "totDist2")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		double parameter = pEffectInfo->internalEffectParameter();
		AlterFunction * pChangeFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, false, true);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, true, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "simDist2")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new CovariateDistance2SimilarityNetworkFunction(networkName,
				covariateName, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2SimilarityNetworkFunction(networkName,
				covariateName, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "simEgoDist2")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new CovariateDistance2EgoAltSimNetworkFunction(networkName,
				covariateName, false, true);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2EgoAltSimNetworkFunction(networkName,
				covariateName, true, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "simEgoInDist2")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new CovariateDistance2EgoAltSimNetworkFunction(networkName,
				covariateName, false, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2EgoAltSimNetworkFunction(networkName,
				covariateName, true, false);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "altInDist2")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new CovariateDistance2InAlterNetworkFunction(networkName,
				covariateName, false, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2InAlterNetworkFunction(networkName,
				covariateName, true, false);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "totInDist2")
	{
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new CovariateDistance2InAlterNetworkFunction(networkName,
				covariateName, false, true);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2InAlterNetworkFunction(networkName,
				covariateName, true, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "altDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		double parameter = pEffectInfo->internalEffectParameter();
		AlterFunction * pChangeFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, false, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, true, false);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "totDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		double parameter = pEffectInfo->internalEffectParameter();
		AlterFunction * pChangeFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, false, true);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, true, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "simDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		AlterFunction * pChangeFunction =
			new CovariateDistance2SimilarityNetworkFunction(networkName,
				covariateName, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2SimilarityNetworkFunction(networkName,
				covariateName, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "altInDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		AlterFunction * pChangeFunction =
			new CovariateDistance2InAlterNetworkFunction(networkName,
				covariateName, false, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2InAlterNetworkFunction(networkName,
				covariateName, true, false);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "totInDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		AlterFunction * pChangeFunction =
			new CovariateDistance2InAlterNetworkFunction(networkName,
				covariateName, false, true);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2InAlterNetworkFunction(networkName,
				covariateName, true, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "simEgoDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		AlterFunction * pChangeFunction =
			new CovariateDistance2EgoAltSimNetworkFunction(networkName,
				covariateName, false, true);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2EgoAltSimNetworkFunction(networkName,
				covariateName, true, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "simEgoInDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		AlterFunction * pChangeFunction =
			new CovariateDistance2EgoAltSimNetworkFunction(networkName,
				covariateName, false, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2EgoAltSimNetworkFunction(networkName,
				covariateName, true, false);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "intercept")
	{
		pEffect = new InterceptEffect(pEffectInfo);

	}
	else if (effectName == "feedback")
	{
		pEffect = new FeedbackEffect(pEffectInfo);
	}
	else if (effectName == "wiener")
	{
		pEffect = new WienerEffect(pEffectInfo);
	}
	else
	{
		throw domain_error("Unexpected effect name: " + effectName);
	}
	return pEffect;
}

}
