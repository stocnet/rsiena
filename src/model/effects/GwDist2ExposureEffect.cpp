/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwDist2ExposureEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwDist2ExposureEffect.
 *****************************************************************************/

#include <stdexcept>
#include <cmath>
#include "GwDist2ExposureEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"

using namespace std;

namespace siena
{

/**
 * Constructor. Reads the internal effect parameter directly from
 * pEffectInfo into its own member (mirroring GwdspExposureEffect and
 * TotalGwInAltDist2NCEffect), rather than from a base-class member,
 * since pEffectInfo is already available here. Validates nonnegativity
 * immediately, same as the analogous behavior effects.
 */
GwDist2ExposureEffect::GwDist2ExposureEffect(const EffectInfo* pEffectInfo) :
    DiffusionRateEffect(pEffectInfo)
{
    this->lparameter = pEffectInfo->internalEffectParameter();
    this->lweight = -0.01 * this->lparameter;
    this->lexpmweight = exp(-this->lweight);
    this->lexpfactor = (1 - exp(this->lweight));
    // Currently not resolved in the effect factory but here at construction time.
    // gwdist2FFExposure (FF, I->K->J): h reached via ties OUT OF j.
    // gwdist2FBExposure (FB, I->K<-J, default): h reached via ties INTO j.
    this->lforward = (this->leffectName == "gwdist2FFExposure");

    if (this->lparameter < 0)
    {
        throw runtime_error(
            "GWDist2 Exposure must have nonnegative internal effect parameter");
    }
}

/**
 * Initializes this effect, then re-validates that the internal effect
 * parameter (the GW decay rate alpha for this effect) is nonnegative,
 * matching the precedent set by the analogous behavior effects
 * (TotalGwdspAlterEffect, TotalGwInAltDist2NCEffect).
 */
void GwDist2ExposureEffect::initialize(const Data* pData,
    State* pState,
    int period,
    Cache* pCache)
{
    DiffusionRateEffect::initialize(pData, pState, period, pCache);

    // The value-sum index (see proximityValue()) is only well defined for
    // non-negative behavior: a negative value would push the accumulated
    // sum below zero and index the lookup out of bounds. Fail loudly here
    // rather than silently produce garbage during simulation.
    if (this->min() < 0)
    {
        throw runtime_error(
            "GWDIST2 exposure requires a non-negative behavior variable");
    }

    // The per-gateway total (see proximityValue()) sums the int-valued
    // behavior scores of the actors reached through gateway j. Those
    // actors are always in the sender node set (that is where the
    // behavior is defined): FB reaches them via inTies(j), FF via
    // outTies(j), and either way there are at most n() of them, each with
    // value <= max(). The largest possible index is therefore n() * max();
    // sizing to n() * max() + 1 makes every reachable index valid, so no
    // run-time clamping is needed. This bound is correct for bipartite
    // (n() != m()) as well as one-mode/symmetric. Contrast GwdspExposure
    // Effect, whose index is a shared-partner COUNT bounded by m().
    double pow_ = 1;
    int size = this->pNetwork()->n() * this->max() + 1;
    this->lcumulativeWeight.resize(size); // default values 0
    for (int i = 1; i < size; i++)
    {
        pow_ *= this->lexpfactor;
        this->lcumulativeWeight[i] = this->lexpmweight * (1 - pow_);
    }
}

/**
 * GWDIST2 (FB and FF): s_i = Sum_j x_ij * gwWeight( Sum_{h!=i} value(h) ),
 * where the inner sum over h is taken via ties INTO j (FB,
 * gwdist2FBExposure, default) or OUT OF j (FF, gwdist2FFExposure).
 *
 * The saturating transform applies once per gateway j, after that
 * gateway's inner value-sum is complete -- one lcumulativeWeight lookup
 * per j, not per h. This is the structurally different kernel from GWDSP
 * (GwdspExposureEffect): GWDIST2 saturates the accumulated value behind
 * a single gateway; GWDSP saturates the path COUNT to a fixed
 * destination, accumulated across all gateways. They cannot share a loop
 * skeleton, which is the other reason these are separate classes -- but
 * both totals are exact non-negative integer indices, so both can use the
 * same precomputed-lookup mechanism (see initialize()).
 *
 * No bounds clamp: initialize() sizes lcumulativeWeight so that every
 * reachable value-sum (0 .. n()*max()) is a valid index, and the
 * non-negativity of the behavior is enforced there too.
 */
double GwDist2ExposureEffect::proximityValue(const Network* pNetwork, int i)
{
    double totalAlterValue = 0;

    if (pNetwork->outDegree(i) > 0)
    {
        for (IncidentTieIterator iter = pNetwork->outTies(i); iter.valid(); iter.next())
        {
            int j = iter.actor();
            int totalAlterInDist2Value = 0;

            IncidentTieIterator iterH = this->lforward ?
                pNetwork->outTies(j) : pNetwork->inTies(j);

            for (; iterH.valid(); iterH.next())
            {
                if (i != iterH.actor())
                {
                    totalAlterInDist2Value += this->value(iterH.actor());
                }
            }

            totalAlterValue += this->lcumulativeWeight[totalAlterInDist2Value];
        }
    }

    return totalAlterValue;
}

}
