/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DoubleCovariateCatFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * DoubleCovariateCatFunction.
 *****************************************************************************/

#include "DoubleCovariateCatFunction.h"
#include <math.h> /* round, sqrt */
#include "utils/Utils.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "network/TieIterator.h"
#include "model/State.h"


using namespace std;

namespace siena
{

/**
 * Creates a new predicate.
 */
DoubleCovariateCatFunction::DoubleCovariateCatFunction(std::string covariateName1,
				std::string covariateName2, 
				std::string networkName, double parameter, bool excludeMissing,
				bool byTies): 
			DoubleCovariateFunction(covariateName1, covariateName2)
{
	this->lpNetwork = 0;
	this->lnetworkName = networkName;
//	this->lpNetworkCache = 0;
	this->lroot = (parameter == 2)||(parameter == 4);
	this->laverage = (parameter >= 3);
	this->lexcludeMissing = excludeMissing;
	this->lbyTies = byTies;

//	this->lFirstCovariateName = covariateName1;
//	this->lSecondCovariateName = covariateName2;
//	this->lpFirstConstantCovariate = 0;
//	this->lpSecondConstantCovariate = 0;
//	this->lpFirstChangingCovariate = 0;
//	this->lpSecondChangingCovariate = 0;	
	
	this->lpTotalCovariateCombinations = 0;
	this->lpNumberTieValues = 0;
	this->lpFirstCovariateNumbers = 0;
	this->lpSecondCovariateNumbers = 0;	
}


/**
 * Initializes this predicate.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DoubleCovariateCatFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	DoubleCovariateFunction::initialize(pData, pState, period, pCache);
	this->lpNetwork = pState->pNetwork(this->lnetworkName); 
// see inline methods in DoubleCovariateCatFunction.h

	int firstMin = int(round(this->firstCovariateMinimum()));
	int firstMax = int(round(this->firstCovariateMaximum()));
	int secondMin = int(round(this->secondCovariateMinimum()));
	this->lSecondMax = int(round(this->secondCovariateMaximum()));
	if (firstMin < 0)
	{
		throw logic_error("sameXsameV: minimum of first covariate is negative");		
	}
	if (secondMin < 0)
	{
		throw logic_error("sameXsameV: minimum of second covariate is negative");		
	}
	if (firstMax > 20)
	{
		throw logic_error("sameXsameV: first covariate has a maximum which is too large");		
	}
	if (this->lSecondMax > 20)
	{
		throw logic_error("sameXsameV: second covariate has a maximum which is too large");		
	}	

	firstMax++;
	this->lSecondMax++;

	this->lpNumberTieValues = new int[this->lSecondMax] {};
	this->lpTotalCovariateCombinations = new int[this->lSecondMax] {};
	this->lpFirstCovariateNumbers = new int[firstMax] {};
	this->lpSecondCovariateNumbers = new int[this->lSecondMax] {};
		
	for (int i = 0; i < this->firstCovariateN(); i++)
	{
		this->lpFirstCovariateNumbers[this->firstCovariateIntValue(i)]++;		
	}
	
	for (int i = 0; i < this->secondCovariateN(); i++)
	{
		this->lpSecondCovariateNumbers[this->secondCovariateIntValue(i)]++;		
	}
}


/**
 * Deallocates this double covariate object.
 */
DoubleCovariateCatFunction::~DoubleCovariateCatFunction()
{
	delete[] this->lpSecondCovariateNumbers;
	this->lpSecondCovariateNumbers = 0;
	delete[] this->lpFirstCovariateNumbers;
	this->lpFirstCovariateNumbers = 0;
	delete[] this->lpTotalCovariateCombinations;
	this->lpTotalCovariateCombinations = 0;
	delete[] this->lpNumberTieValues;
	this->lpNumberTieValues = 0;
}

void DoubleCovariateCatFunction::preprocessEgo(int ego)
{
	DoubleCovariateFunction::preprocessEgo(ego);
	
	for (int b = 0; b < this->lSecondMax; b++)
	{
		this->lpNumberTieValues[b]= 0;
	}
	
	if ((!this->firstMissing(ego))||(!this->lexcludeMissing))
	{	
		if (this->lbyTies) // iterate by ties
		{
			if (this->lexcludeMissing)
			{			
				for (TieIterator iter = this->pNetwork()->ties();
										iter.valid(); iter.next())
				{
					// Get the sender and receiver of the  tie.
					int i = iter.ego();
					int h = iter.alter();
				// i needs to have the same covariate value as ego:
					if (!(this->firstMissing(i)))
					{
						if ((this->firstCovariateIntValue(i) == this->firstCovariateIntValue(ego)) &&
											(!this->secondMissing(h)))
						{
							this->lpNumberTieValues[this->secondCovariateIntValue(h)]++;
						}
					}
				}
			}
			else
			{		
				for (TieIterator iter = this->pNetwork()->ties();
										iter.valid(); iter.next())
				{
				// iter.ego needs to have the same covariate value as ego:
					if (this->firstCovariateIntValue(iter.ego()) == this->firstCovariateIntValue(ego))
					{
						this->lpNumberTieValues[this->secondCovariateIntValue(iter.alter())]++;
					}
				}
			}
		}
		else // iterate by nodes and outTies
		{
			if (this->lexcludeMissing)
			{
				for (int i = 0; i < this->firstCovariateN(); i++)
				{
				// i needs to have the same covariate value as ego:
					if (!(this->firstMissing(i) || this->firstMissing(ego)))
					{
						if (this->firstCovariateIntValue(i) == this->firstCovariateIntValue(ego))
						{
							for (IncidentTieIterator iter =	this->pNetwork()->outTies(i);
							iter.valid();
							iter.next())
							{
						// Get the receiver of the outgoing tie.
								int h = iter.actor();
								if (!this->secondMissing(h))
								{	
									this->lpNumberTieValues[this->secondCovariateIntValue(h)]++;
								}
							}
						}
					}
				}
			}
			else
			{		
				for (int i = 0; i < this->firstCovariateN(); i++)
				{
					// i needs to have the same covariate value as ego:	
					if (this->firstCovariateIntValue(i) == this->firstCovariateIntValue(ego))
					{
						for (IncidentTieIterator iter =	this->pNetwork()->outTies(i);
						iter.valid();
						iter.next())
						{
							this->lpNumberTieValues[this->secondCovariateIntValue(iter.actor())]++;
						}
					}	
				}
			}
		}
	}
}

/**
 * Returns the numbers of nodes with value a for the first covariate.
 */
int DoubleCovariateCatFunction::firstCovariateNumbers(int a) const
{
	return lpFirstCovariateNumbers[a];
}

/**
 * Returns the numbers of nodes with value b for the second covariate.
 */
int DoubleCovariateCatFunction::secondCovariateNumbers(int b) const
{
	return lpSecondCovariateNumbers[b];
}

/**
 * Returns the value of the covariate for the given pair of actors.
 */
int DoubleCovariateCatFunction::numberCovariateTies(int b) const
{
	return this->lpNumberTieValues[b];
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double DoubleCovariateCatFunction::value(int alter) const
{
	double value = 0;
	if ((this->firstCovariateIntValue(this->ego()) >= 1) &&
			(this->secondCovariateIntValue(alter) >= 1))
	{
		value = numberCovariateTies(this->secondCovariateIntValue(alter));
		if (this->laverage)
		{
			value /= (this->firstCovariateNumbers(this->firstCovariateIntValue(this->ego())) *
					this->secondCovariateNumbers(this->secondCovariateIntValue(alter)));
// the denominator cannot be 0
		}
		if (this->lroot)
		{
			value = sqrt(value);
		}
	}
	return value;
}

}
