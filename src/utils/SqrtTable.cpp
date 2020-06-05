/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SqrtTable.cpp
 *
 * Description: This file contains the implementation of the
 * SqrtTable class.
 *****************************************************************************/

#include <cmath>

#include "SqrtTable.h"

namespace siena
{

// The size of the lookup table. The square roots of integers larger than
// this will be calculated with sqrt() directly.

const int LIMIT = 1000;


// Nullify the pointer to the single instance so we can allocate it
// on demand.

SqrtTable * SqrtTable::lpInstance = 0;


/**
 * Returns the only instance of the SqrtTable class. The instance
 * is allocated on demand.
 */
SqrtTable * SqrtTable::instance()
{
	if (!SqrtTable::lpInstance)
	{
		SqrtTable::lpInstance = new SqrtTable();
	}

	return SqrtTable::lpInstance;
}


/**
 * Constructor.
 */
SqrtTable::SqrtTable()
{
	this->ltable = new double[LIMIT];

	// Leave the entries undefined, so they can be computed on demand.

	for (int i = 0; i < LIMIT; i++)
	{
		this->ltable[i] = -1;
	}
}


/**
 * Destructor.
 */
SqrtTable::~SqrtTable()
{
	delete[] this->ltable;
}


/**
 * Returns the square root of the given integer.
 */
double SqrtTable::sqrt(int i)
{
	double root = 0;

	if (i < LIMIT)
	{
		if (this->ltable[i] < 0)
		{
			this->ltable[i] = std::sqrt((double) i);
		}

		root = this->ltable[i];
	}
	else
	{
		root = std::sqrt((double) i);
	}

	return root;
}

}
