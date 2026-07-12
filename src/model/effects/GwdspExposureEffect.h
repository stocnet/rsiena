/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwdspExposureEffect.h
 *
 * Description: This file contains the definition of the
 * GwdspExposureEffect class, which handles the geometrically
 * weighted dyadwise shared partners (GWDSP) diffusion rate effects
 * (gwdspFBExposure, gwdspFFExposure).
 *
 * Split out from Distance2ExposureEffect: GWDSP saturates the COUNT of
 * two-paths between ego and a fixed distance-2 alter h (accumulated
 * across all gateways j before any saturation is applied), which is a
 * structurally different traversal from the GW kernel that
 * Distance2ExposureEffect's other variants share (saturate a value-sum
 * behind a single gateway j) -- see GwDist2ExposureEffect. Both totals
 * are exact non-negative integer counts, so both classes independently
 * use the same precomputed lcumulativeWeight lookup mechanism, each over
 * its own domain (path count here vs. gateway value-sum there); there is
 * no shared gwWeight() on the base class.
 *****************************************************************************/

#ifndef GWDSPEXPOSUREEFFECT_H_
#define GWDSPEXPOSUREEFFECT_H_

#include "DiffusionRateEffect.h"
#include <vector>

namespace siena
{

class GwdspExposureEffect : public DiffusionRateEffect
{
public:
    GwdspExposureEffect(const EffectInfo* pEffectInfo);

    virtual void initialize(const Data* pData, State* pState, int period, Cache* pCache) override;


protected:
    virtual double proximityValue(const Network* pNetwork, int i) override;

private:
	std::vector<double> lcumulativeWeight;
	double lparameter {};
	double lforward {};
	double lweight {};
	double lexpmweight {};
	double lexpfactor {};

	std::vector<int> lTwoPathCount;
	std::vector<int> lTouched;

};

}

#endif /* GWDSPEXPOSUREEFFECT_H_ */
