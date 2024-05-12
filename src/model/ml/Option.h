#ifndef OPTION_H_
#define OPTION_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * This class represents the option of a ministep.
 */
class Option
{
public:
	Option(int variableIndex, int ego, int alter = 0);

	inline int variableIndex() const;
	inline int ego() const;
	inline int alter() const;

private:
	// The index of the respective variable
	int lvariableIndex{};

	// The ego
	int lego{};

	// The alter (0 for behavior variable)
	int lalter{};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the index of the variable corresponding to this option.
 */
int Option::variableIndex() const
{
	return this->lvariableIndex;
}


/**
 * Returns the ego of this option.
 */
int Option::ego() const
{
	return this->lego;
}


/**
 * Returns the alter of this option (0, if the respective variable
 * is a behavior variable).
 */
int Option::alter() const
{
	return this->lalter;
}


// ----------------------------------------------------------------------------
// Section: Operators
// ----------------------------------------------------------------------------

bool operator<(const Option & rOption1, const Option & rOption2);

bool operator==(const Option & rOption1, const Option & rOption2);
}

#endif /* OPTION_H_ */
