/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * TransitiveTripletsEffect.
 *****************************************************************************/

#ifndef AGREEMENT_TRANSITIVITY_EFFECT_H_
#define AGREEMENT_TRANSITIVITY_EFFECT_H_

#include "network/Network.h"
#include "model/effects/NetworkEffect.h"

namespace siena {

/**
 * Real transitivity effect.
 *
 * <pre>
 *       h                         h
 *      ^ ^                       ^ ^
 *     /   \          ------->   /   \
 *  ego     alter             ego --> alter
 * </pre>
 *
 * The agreement (ego, alter -> h) existed, the transitive edge (ego -> alter)
 * flipped into existence.
 */
class AgreementTransitivityEffect: public NetworkEffect {
public:
	AgreementTransitivityEffect(const EffectInfo * pEffectInfo) :
			NetworkEffect(pEffectInfo) {
	}
	double calculateContribution(int alter) const {
		// GMM effect, don't need contribution implementation
		return 0;
	}
protected:
	double tieStatistic(int alter) {
		int statistic = 0;
		// observation form the start of the period
		const Network* pStart = this->pData()->pNetwork(this->period());
		// (ego -> alter) exists in the final state (call to tieStatistics)
		// (ego -> alter) missing in the start state
		if (pStart->tieValue(ego(), alter) == 0) {
			// loop over all (ego -> h), in the start state
			for (IncidentTieIterator h = pStart->outTies(ego()); h.valid();
					h.next()) {
				// check existence in the final state
				if (outTieExists(h.actor())) {
					// count the (h -> alter) in both states
					statistic += pNetwork()->tieValue(alter, h.actor())
							* pStart->tieValue(alter, h.actor());
				}
			}
		}
		return statistic;
	}
};

}

#endif // AGREEMENT_TRANSITIVITY_EFFECT_H_
