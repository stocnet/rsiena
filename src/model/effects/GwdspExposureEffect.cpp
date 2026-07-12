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

    this->lTwoPathCount.assign(this->pNetwork()->n(), 0);
    this->lTouched.clear();
}

/**
 * GWDSP (FB and FF): s_i = Sum_{h!=i} value(h) * lcumulativeWeight[pathCount(i,h)].
 *
 * pathCount(i,h) = Sum_j x_ij*x_hj (FB, gwdspFBExposure) or
 * Sum_j x_ij*x_jh (FF, gwdspFFExposure) -- the number of two-paths
 * from i to h.
 *
 * Rather than looping h over all n actors and asking pNetwork->twoPathCount/
 * inTwoStarCount(i, h) from scratch for each (an O(n) sweep, each query an
 * O(deg) merge-walk -- expensive because it interrogates every possible h,
 * including the vast majority sharing nothing with i), this accumulates the
 * same counts by summing over i's actual out-ties instead:
 *
 *   pathCount(i,h) = Sum_j [i->j] * [h->j]           (FB)
 *                  = Sum_j [i->j] * [j->h]           (FF)
 */
double GwdspExposureEffect::proximityValue(const Network* pNetwork, int i)
{
    for (IncidentTieIterator iter = pNetwork->outTies(i); iter.valid(); iter.next())
    {
        int j = iter.actor();

        IncidentTieIterator iterH = this->lforward ?
            pNetwork->outTies(j) : pNetwork->inTies(j);

        for (; iterH.valid(); iterH.next())
        {
            int h = iterH.actor();
            if (h == i)
            {
                continue;
            }
            if (this->lTwoPathCount[h] == 0)
            {
                this->lTouched.push_back(h);
            }
            this->lTwoPathCount[h]++;
        }
    }

    double totalAlterValue = 0;
    // No bounds clamp: initialize() sizes lcumulativeWeight to m() + 1,
    // so every reachable partner count (0 .. m()) is a valid index.
    for (int h : this->lTouched)
    {
        totalAlterValue += this->value(h) * this->lcumulativeWeight[this->lTwoPathCount[h]];
        this->lTwoPathCount[h] = 0;
    }
    this->lTouched.clear();

    return totalAlterValue;
}

}
