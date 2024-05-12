#ifndef COVARIATEPREDICATE_H_
#define COVARIATEPREDICATE_H_

#include "AlterPredicate.h"
#include "utils/NamedObject.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;
class ChangingCovariate;
class BehaviorLongitudinalData;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines an alter predicate that depends on a certain covariate or
 * behavior variable.
 */
class CovariatePredicate: public AlterPredicate, NamedObject {
public:
	CovariatePredicate(std::string covariateName);
	virtual ~CovariatePredicate() {}
	virtual void initialize(const Data * pData, State * pState, int period,
			Cache * pCache);

protected:
	double covariateValue(int i) const;
	bool missing(int i) const;

private:
	ConstantCovariate * lpConstantCovariate;
	ChangingCovariate * lpChangingCovariate;
	BehaviorLongitudinalData * lpBehaviorData;

	// The current value of a behavior variable per each actor.
	// This array is 0 for covariate-based effects.

	const int * lvalues {};
};

}

#endif /* COVARIATEPREDICATE_H_ */
