/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DiffusionRateEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * DiffusionRateEffect.
 *****************************************************************************/

#include <cmath>
#include <cstring>
#include "DiffusionRateEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/IncidentTieIterator.h"
#include <R_ext/Print.h>

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pNetwork the network this effect depends on
 * @param[in] pBehaviorVariable the behavior variable this effect depends on
 * @param[in] effectName the name of this effect
 * @param[in] parameter the statistical parameter of this effect
 * @param[in] internalEffectParameter the internal effect parameter
 */
DiffusionRateEffect::DiffusionRateEffect(const EffectInfo* pEffectInfo,
    const Network * pNetwork,
    const BehaviorVariable * pBehaviorVariable,
    string effectName,
    double parameter,
    double internalEffectParameter) : Effect(pEffectInfo)
{
    this->lpNetwork = pNetwork;
    this->lpBehaviorVariable = pBehaviorVariable;
    this->lpConstantCovariate = 0;
    this->lpChangingCovariate = 0;
    this->leffectName = effectName;
    this->lparameter = parameter;
    this->linternalEffectParameter = round(internalEffectParameter);
    this->labsInternalEffectParameter = std::abs(this->linternalEffectParameter);
    this->linternalNonZero = (this->linternalEffectParameter != 0);

    if (((effectName == "infectDeg") || (effectName == "infectIn") ||
                (effectName == "infectOut")) && (this->linternalEffectParameter < 0))
    {
        throw logic_error("Negative internal parameter not permitted for effect "+effectName);
    }
}

/**
 * Constructor.
 * @param[in] pNetwork the network this effect depends on
 * @param[in] pBehaviorVariable the behavior variable this effect depends on
 * @param[in] pConstantCovariate the covariate this effect depends on
 * @param[in] pChangingCovariate the changing covariate this effect depends on
 * @param[in] effectName the name of this effect
 * @param[in] parameter the statistical parameter of this effect
 * @param[in] internalEffectParameter the internal effect parameter
 */
DiffusionRateEffect::DiffusionRateEffect(const EffectInfo* pEffectInfo,
    const Network * pNetwork,
    const BehaviorVariable * pBehaviorVariable,
    const ConstantCovariate * pConstantCovariate,
    const ChangingCovariate * pChangingCovariate,
    string effectName,
    double parameter,
    double internalEffectParameter) : Effect(pEffectInfo)
{
    this->lpNetwork = pNetwork;
    this->lpBehaviorVariable = pBehaviorVariable;
    this->lpChangingCovariate = pChangingCovariate;
    this->lpConstantCovariate = pConstantCovariate;
    this->leffectName = effectName;
    this->lparameter = parameter;
    this->linternalEffectParameter = round(internalEffectParameter);
    this->labsInternalEffectParameter = std::abs(this->linternalEffectParameter);
    this->linternalNonZero = (this->linternalEffectParameter != 0);

    if ((effectName == "infectCovar") && (this->linternalEffectParameter < 0))
    {
        throw logic_error("Negative internal parameter not permitted for effect "+effectName);
    }
}

/**
 * Destructor.
 */
DiffusionRateEffect::~DiffusionRateEffect()
{
    // Nothing to clean up
}

double DiffusionRateEffect::proximityValue(const Network* pNetwork, int i, int period) const
{
    throw std::logic_error("proximityValue not implemented for this effect type.");
}
/**
 * Returns the raw social proximity/exposure statistic for actor i.
 * This calculates the statistic according to the formula:
 * a_i(y) = sum_j (z_j * d_ij(y))
 * where d_ij(y) is a measure of closeness/influence between actors i and j,
 * and z_j indicates whether j has adopted the innovation.
 * 
 * ALL conditional logic for different effect types is contained here.
 */
/* double DiffusionRateEffect::proximityValue(const Network * pNetwork, int i,
        int period) const
{
    int egoNumer = 1;
    int egoDenom = 1;

    // Determine numerator and denominator based on effect type
    if (this->leffectName == "avExposure" || 
        this->leffectName == "avTinExposureDist2")
    {
        egoDenom = max(1, pNetwork->outDegree(i));
    }
    else if (this->leffectName == "susceptAvIn")
    {
        egoNumer = pNetwork->inDegree(i);
        egoDenom = max(1, pNetwork->outDegree(i));
    }

    // Calculate the sum of proximity measures
    double totalAlterValue = 0;
    int numInfectedAlter = 0;
    
    if (pNetwork->outDegree(i) > 0)
    {
        for (IncidentTieIterator iter = pNetwork->outTies(i);
             iter.valid();
             iter.next())
        {
            // Distance-2 effects: sum over i's alters' in-neighbors
            if (leffectName == "anyInExposureDist2" || 
                leffectName == "totInExposureDist2" || 
                leffectName == "avTinExposureDist2" || 
                leffectName == "totAInExposureDist2")
            {
                int j = iter.actor();
                double totalAlterInDist2Value = 0;
                
                for (IncidentTieIterator iterH = pNetwork->inTies(j);
                    iterH.valid();
                    iterH.next())
                {	
                    if (i != iterH.actor())
                    {
                        double alterInDist2Value = lpBehaviorVariable->value(iterH.actor());
                        if (alterInDist2Value >= 0.5)
                        {
                            numInfectedAlter++;
                        }
                        totalAlterInDist2Value += alterInDist2Value;
                    }
                }
                
                // Apply distance-2 specific normalizations
                if ((leffectName == "totAInExposureDist2") && 
                    ((pNetwork->inDegree(j) - 1) > 0))
                {
                    totalAlterInDist2Value /= (pNetwork->inDegree(j) - 1);
                }
                if (leffectName == "anyInExposureDist2")
                {
                    totalAlterInDist2Value = std::min(totalAlterInDist2Value, 1.0);
                }
                
                totalAlterValue += totalAlterInDist2Value;
            }
            // Standard proximity effects (distance-1)
            else if (this->leffectName == "totExposure" ||
                     this->leffectName == "avExposure" ||
                     this->leffectName == "susceptAvIn" ||
                     this->leffectName == "infectDeg" ||
                     this->leffectName == "infectIn" ||
                     this->leffectName == "infectOut" ||
                     this->leffectName == "susceptAvCovar")
            {
                double alterValue = this->lpBehaviorVariable->value(iter.actor());
                
                if (alterValue >= 0.5)
                {
                    numInfectedAlter++;
                }
                
                // Weighted by alter's degree
                if (this->leffectName == "infectIn")
                {
                    alterValue *= pNetwork->inDegree(iter.actor());
                }
                else if ((this->leffectName == "infectDeg") ||
                         (this->leffectName == "infectOut"))
                {
                    alterValue *= pNetwork->outDegree(iter.actor());
                }
                
                totalAlterValue += alterValue;
            }
            // Covariate-weighted proximity
            else if (this->leffectName == "infectCovar")
            {
                double alterValue = this->lpBehaviorVariable->value(iter.actor());

                if (this->lpConstantCovariate)
                {
                    alterValue *= this->lpConstantCovariate->value(iter.actor());
                }
                else if (this->lpChangingCovariate)
                {
                    alterValue *= this->lpChangingCovariate->value(iter.actor(), period);
                }
                else
                {
                    throw logic_error("No individual covariate found for infectCovar.");
                }

                totalAlterValue += alterValue;
            }
        }
    }
    
    // Apply internal effect parameter thresholding
    if (this->linternalNonZero)
    {
        if (numInfectedAlter < this->labsInternalEffectParameter)
        {
            totalAlterValue = 0;
        }
        else if (this->linternalEffectParameter < 0)
        {
            if (totalAlterValue > this->labsInternalEffectParameter)
            {
                totalAlterValue = this->labsInternalEffectParameter;
            }
        }
    }

    // Apply ego susceptibility weighting
    totalAlterValue *= egoNumer;
    
    // Normalize by denominator
    double rawStatistic;
    if (egoDenom > 1)
    {
        rawStatistic = totalAlterValue / egoDenom;
    }
    else
    {
        rawStatistic = totalAlterValue;
    }
    
    // Apply ego covariate susceptibility weighting
    if (this->leffectName == "susceptAvCovar")
    {
        if (this->lpConstantCovariate)
        {
            rawStatistic *= this->lpConstantCovariate->value(i);
        }
        else if (this->lpChangingCovariate)
        {
            rawStatistic *= this->lpChangingCovariate->value(i, period);
        }
        else
        {
            throw logic_error("No individual covariate found for susceptAvCovar.");
        }
    }
    
    return rawStatistic;
} */

double DiffusionRateEffect::applyInternalEffectParameter(double value, int numInfectedAlter) const
{
    if (this->linternalNonZero)
    {
        if (numInfectedAlter < this->labsInternalEffectParameter)
        {
            value = 0;
        }
        else if (this->linternalEffectParameter < 0)
        {
            if (value > this->labsInternalEffectParameter)
            {
                value = this->labsInternalEffectParameter;
            }
        }
    }
    return value;
}

/**
 * Returns the raw statistic (for scores and statistics calculation).
 * Simple transformer: just returns the proximity value.
 */
double DiffusionRateEffect::value(int i, int period) const
{
    return this->proximityValue(this->lpNetwork, i, period);
}

/**
 * Returns the exponentiated rate contribution (for rate calculations).
 * Exponentiates the raw statistic: exp(parameter * rawStatistic)
 */
double DiffusionRateEffect::logRate(int i, int period) const
{
    double rawStatistic = this->proximityValue(this->lpNetwork, i, period);    
    return this->lparameter * rawStatistic;
}

/**
 * Stores the parameter for the diffusion rate effect.
 */
void DiffusionRateEffect::parameter(double parameterValue)
{
    this->lparameter = parameterValue;
}

/**
 * Returns the parameter for the diffusion rate effect.
 */
double DiffusionRateEffect::parameter() const
{
    return this->lparameter;
}

/**
 * Stores the internal effect parameter for the diffusion rate effect.
 */
void DiffusionRateEffect::setInternalEffectParameter(int parValue)
{
    this->linternalEffectParameter = parValue;
    this->labsInternalEffectParameter = std::abs(this->linternalEffectParameter);
    this->linternalNonZero = (this->linternalEffectParameter != 0);
}

/**
 * Returns the internal effect parameter for the diffusion rate effect.
 */
int DiffusionRateEffect::getInternalEffectParameter() const
{
    return this->linternalEffectParameter;
}

}