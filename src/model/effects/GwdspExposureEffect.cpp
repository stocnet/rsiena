/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwdspExposureEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwdspExposureEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "GwdspExposureEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
 
using namespace std;
 
namespace siena
{
 
/**
 * Constructor. Computes the GW decay parameters directly from pEffectInfo.
 * Validates nonnegativity immediately, same as the
 * analogous behavior effect.
 */
GwdspExposureEffect::GwdspExposureEffect(const EffectInfo* pEffectInfo) :
    DiffusionRateEffect(pEffectInfo)
{
    this->lparameter  = pEffectInfo->internalEffectParameter();
    this->lweight = -0.01 * lparameter;
    this->lexpmweight = exp(-this->lweight);
    this->lexpfactor = (1 - exp(this->lweight));
    this->lforward = (this->leffectName == "gwdspFFExposure"); 
    
    if (lparameter < 0)
    {
        throw runtime_error(
            "GWDSP Exposure must have nonnegative internal effect parameter");
    }
}

void GwdspExposureEffect::initialize(const Data* pData,
    State* pState,
    int period,
    Cache* pCache)
{
    DiffusionRateEffect::initialize(pData, pState, period, pCache);

    // The index is a shared-partner count between two senders: the number
    // of receivers both point to (inTwoStarCount, FB) or the number of
    // two-paths (twoPathCount, FF). Either is bounded by the receiver set,
    // so the largest possible index is m() (two senders sharing every
    // receiver -- reachable on a bipartite network). Sizing to m() + 1
    // makes every reachable index valid, so no run-time clamping is needed.
    double pow_ = 1;
    int size = this->pNetwork()->m() + 1;
    this->lcumulativeWeight.resize(size); // default values 0
    for (int i = 1; i < size; i++)
    {
        pow_ *= this->lexpfactor;
        this->lcumulativeWeight[i] = this->lexpmweight * (1 - pow_);
    }
}

/**
 * GWDSP (FB and FF): s_i = Sum_{h!=i} value(h) * lcumulativeWeight[pathCount(i,h)].
 *
 * pathCount(i,h) = Sum_j x_ij*x_hj (FB, gwdspFBExposure) or
 * Sum_j x_ij*x_jh (FF, gwdspFFExposure) -- the number of two-paths
 * from i to h. Unlike the GWDIST2 kernel, the
 * saturating transform here applies per DESTINATION h, accumulated across
 * every gateway j, so it cannot be computed inside a single pass over j;
 * h's count must be complete before the lookup is applied to it.
 */
double GwdspExposureEffect::proximityValue(const Network* pNetwork, int i) const
{

    double totalAlterValue = 0;
    int n = pNetwork->n();

    // Loop over all N actors in the network
    for (int h = 0; h < n; h++)
    {
        if (h == i) 
        {
            continue;
        }

        int pathCount = this->lforward ? pNetwork->twoPathCount(i, h)
                                       : pNetwork->inTwoStarCount(i, h);

        // No bounds clamp: initialize() sizes lcumulativeWeight to m() + 1,
        // so every reachable partner count (0 .. m()) is a valid index.
        totalAlterValue += this->value(h) * this->lcumulativeWeight[pathCount];
    }

    return totalAlterValue;
}

}
