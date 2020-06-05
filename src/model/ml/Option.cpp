#include "Option.h"

namespace siena
{

/**
 * Constructs a new option with the given properties.
 */
Option::Option(int variableIndex, int ego, int alter)
{
	this->lvariableIndex = variableIndex;
	this->lego = ego;
	this->lalter = alter;
}


/**
 * Returns if this option is less than the given option.
 */
bool operator<(const Option & rOption1, const Option & rOption2)
{
	return rOption1.variableIndex() < rOption2.variableIndex() ||
		(rOption1.variableIndex() == rOption2.variableIndex() &&
			(rOption1.ego() < rOption2.ego() ||
				(rOption1.ego() == rOption2.ego() &&
					rOption1.alter() < rOption2.alter())));
}

/**
 * Returns if this option is equal to the given option.
 */
bool operator==(const Option & rOption1, const Option & rOption2)
{
	return (rOption1.variableIndex() == rOption2.variableIndex() &&
		(rOption1.ego() == rOption2.ego() &&
			rOption1.alter() == rOption2.alter()));
}
}
