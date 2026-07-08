/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwDist2ExposureEffect.h
 *
 * Description: This file contains the definition of the
 * GwDist2ExposureEffect class, which handles the geometrically weighted
 * distance-2 (GWDIST2) diffusion rate effects (gwdist2FBExposure,
 * gwdist2FFExposure).
 *
 * Split out from Distance2ExposureEffect: that class uses
 * internalEffectParameter() as an integer threshold/clip count via
 * applyThreshold(). This class uses the same EffectInfo slot for an
 * unrelated purpose -- the continuous GW decay rate alpha, applied via a
 * precomputed lcumulativeWeight lookup (mirroring GwdspExposureEffect).
 * Note the two classes index that lookup on DIFFERENT quantities: GWDSP
 * indexes on a shared-partner COUNT (bounded by the receiver set, m()),
 * whereas GWDIST2 indexes on a per-gateway VALUE-SUM of the behaviour
 * scores of the actors reached through the gateway -- at most n() such
 * actors, each with value <= max(), so bounded by n() * max(). These
 * coincide only when n() == m() (one-mode/symmetric); the bipartite case
 * disambiguates them. Keeping GWDIST2 and the threshold effects in one
 * class would mean the same parameter silently serving two incompatible
 * readings on every call; this class therefore does not call
 * applyThreshold() at all.
 *****************************************************************************/

#ifndef GWDIST2EXPOSUREEFFECT_H_
#define GWDIST2EXPOSUREEFFECT_H_

#include "DiffusionRateEffect.h"
#include <vector>

namespace siena
{

class GwDist2ExposureEffect : public DiffusionRateEffect
{
public:
    GwDist2ExposureEffect(const EffectInfo* pEffectInfo);

    virtual void initialize(const Data* pData, State* pState, int period, Cache* pCache) override;

protected:
    virtual double proximityValue(const Network* pNetwork, int i) const override;

private:
	std::vector<double> lcumulativeWeight;
    bool lforward {};
    double lparameter {};
  	double lweight {};
	double lexpmweight {};
	double lexpfactor {};
};

}

#endif /* GWDIST2EXPOSUREEFFECT_H_ */
