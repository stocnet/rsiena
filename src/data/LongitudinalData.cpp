/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: LongitudinalData.cpp
 *
 * Description: This file contains the implementation of the
 * LongitudinalData class.
 *****************************************************************************/

#include "LongitudinalData.h"
#include "data/ActorSet.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction and destruction
// ----------------------------------------------------------------------------

/**
 * Creates a new longitudinal data object for the given set of actors and
 * number of observations.
 * @param[in] id the ID that is unique among all longitudinal data object
 * of the parent Data instance
 * @param[in] name the name of the corresponding dependent variable
 */
LongitudinalData::LongitudinalData(int id,
	std::string name,
	const ActorSet * pActorSet,
	int observationCount) : NamedObject(name)
{
	this->lid = id;
	this->lpActorSet = pActorSet;
	this->lobservationCount = observationCount;

	this->lupOnly = new bool[observationCount - 1];
	this->ldownOnly = new bool[observationCount - 1];

	for (int i = 0; i < observationCount - 1; i++)
	{
		this->lupOnly[i] = false;
		this->ldownOnly[i] = false;
	}
}


/**
 * Deallocates this data object.
 */
LongitudinalData::~LongitudinalData()
{
	delete[] this->lupOnly;
	delete[] this->ldownOnly;
	this->lupOnly = 0;
	this->ldownOnly = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the set of actors underlying this data object.
 */
const ActorSet * LongitudinalData::pActorSet() const
{
	return this->lpActorSet;
}


/**
 * Returns the number of observations to be stored in this data object.
 */
int LongitudinalData::observationCount() const
{
	return this->lobservationCount;
}


/**
 * Returns the number of actors for this data object.
 */
int LongitudinalData::n() const
{
	return this->lpActorSet->n();
}


/**
 * Stores if only upward changes were observed at the given period.
 */
void LongitudinalData::upOnly(int period, bool flag)
{
	this->lupOnly[period] = flag;
}


/**
 * Returns if only upward changes were observed at the given period.
 */
bool LongitudinalData::upOnly(int period) const
{
	return this->lupOnly[period];
}


/**
 * Stores if only downward changes were observed at the given period.
 */
void LongitudinalData::downOnly(int period, bool flag)
{
	this->ldownOnly[period] = flag;
}


/**
 * Returns if only downward changes were observed at the given period.
 */
bool LongitudinalData::downOnly(int period) const
{
	return this->ldownOnly[period];
}

}
