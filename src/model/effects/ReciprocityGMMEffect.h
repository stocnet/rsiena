/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ReciprocityGMMEffect.h
 *****************************************************************************/

#ifndef NEW_RECIPROCITY_EFFECT_H_
#define NEW_RECIPROCITY_EFFECT_H_

#include "network/Network.h"
#include "model/effects/NetworkEffect.h"

namespace siena {

/**
 * This class defines the <i>new/real/persistent reciprocity</i> effects, which
 * form the <i>reciprocity</i> family for the GMM model.
 */
class ReciprocityGMMEffect: public NetworkEffect {
public:

	const static int NEW = 0;
	const static int REAL = 1;
	const static int PERSISTENT = 2;

	/**
	 * @param startEdgeSum number of edges need to be present in the start
	 * observation for the dyad to be counted. 0: new, 1: real, 2: persistent reciprocity.
	 */
	ReciprocityGMMEffect(const EffectInfo * pEffectInfo, const int startEdgeSum) :
			NetworkEffect(pEffectInfo) {
		this->startEdgeSum = startEdgeSum;
	}

	double calculateContribution(int) const {
		return 0;
	}

protected:
	virtual double tieStatistic(int alter) {
		// Outgoing tie exists, otherwise we would not be here.
		// Check if it is reciprocated.
		if (this->inTieExists(alter)) {
			// We have reciprocity in the end state.
			// Now look at the start of the period.
			const Network* pStart = this->pData()->pNetwork(this->period());
			if (pStart->tieValue(this->ego(), alter)
					+ pStart->tieValue(alter, this->ego()) == startEdgeSum) {
				return 1;
			}
		}
		return 0;
	}

private:
	int startEdgeSum {};

};

} // namespace siena
#endif // NEW_RECIPROCITY_EFFECT_H_
