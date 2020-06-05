/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DiffusionEffectValueTable.cpp
 *
 * Description: This file contains the implementation of the
 * DiffusionEffectValueTable class.
 *****************************************************************************/
#include <cmath>
#include "DiffusionEffectValueTable.h"

namespace siena {

/**
 * Creates a new look-up table for <i>n</i> actors.
 */

DiffusionEffectValueTable::DiffusionEffectValueTable(int numeratorRange,
		int denominatorRange) {
	this->lvalues = new double[numeratorRange * denominatorRange];
	this->lparameterValues = new double[numeratorRange * denominatorRange];

	// Calculate the table for parameter 0

	this->lparameter = 0;
	this->ldenominatorRange = denominatorRange;
	this->lnumeratorRange = numeratorRange;

	for (int i = 0; i < numeratorRange * denominatorRange; i++) {
		// exp(0) = 1

		this->lvalues[i] = 1;
		this->lparameterValues[i] = 0;
	}
}

/**
 * Deallocates this look-up table.
 */
DiffusionEffectValueTable::~DiffusionEffectValueTable() {
	delete[] this->lvalues;
	delete[] this->lparameterValues;

	this->lvalues = 0;
	this->lparameterValues = 0;
}

/**
 * Returns the current value of the effect parameter.
 */
double DiffusionEffectValueTable::parameter() const {
	return this->lparameter;
}

/**
 * Stores the effect parameter.
 */
void DiffusionEffectValueTable::parameter(double value) {
	this->lparameter = value;
}

/**
 * Returns the value of the effect.
 */

double DiffusionEffectValueTable::value(int numerator, int denominator) {
	int arrayIndex = ((numerator - 1) * this->ldenominatorRange)
			+ (denominator - 1);

	if (this->lparameterValues[arrayIndex] != this->lparameter) {
		// The value stored in the table was calculated for a different
		// parameter, hence we must recalculate the value.

		this->lvalues[arrayIndex] = std::exp(
				this->lparameter * (double) numerator / (double) denominator);
		this->lparameterValues[arrayIndex] = this->lparameter;
	}

	return this->lvalues[arrayIndex];
}

}
