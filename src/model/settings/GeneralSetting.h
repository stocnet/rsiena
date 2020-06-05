/*
 * GeneralSetting.h
 *
 *  Created on: 24.07.2014
 *      Author: ortmann
 */

#ifndef GENERALSETTING_H_
#define GENERALSETTING_H_

#include "Setting.h"

namespace siena {

class GeneralSetting: public Setting {

public:

	virtual ~GeneralSetting();

	ITieIterator* getPermittedSteps();

	int getPermittedSize();

	void initPermittedSteps(const bool* const permitted);

protected:

	GeneralSetting();

	virtual void terminateSetting();

private:

	ITieIterator* lppermittedIter;

};

} /* namespace siena */
#endif /* GENERALSETTING_H_ */
