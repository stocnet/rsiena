/*
 * GeneralSetting.cpp
 *
 *  Created on: 24.07.2014
 *      Author: ortmann
 */

#include <stdexcept>
#include "GeneralSetting.h"

#include "../../network/iterators/UnionTieIterator.h"
#include "../../network/iterators/SingleIterator.h"
#include "../../network/iterators/FilteredIterator.h"


using namespace std;

namespace siena {

GeneralSetting::~GeneralSetting() {
	if (lppermittedIter != 0) {
		delete lppermittedIter;
	}
}

ITieIterator* GeneralSetting::getPermittedSteps() {
	return lppermittedIter->clone();
}

int GeneralSetting::getPermittedSize() {
	return lppermittedIter->size();
}

GeneralSetting::GeneralSetting() :
		Setting(), //
		lppermittedIter(0) {
}

void GeneralSetting::terminateSetting() {
	if (lppermittedIter != 0) {
		delete lppermittedIter;
		lppermittedIter = 0;
	} 
	else 
	{
		throw runtime_error("setting has not been initialized");
	}
}

void GeneralSetting::initPermittedSteps(const bool* const permitted) {
	if (lppermittedIter == 0) {
	ITieIterator* iter = getSteps();
	lppermittedIter = new FilteredIterator(*iter,permitted);
	delete iter;
	} else 
	{
		throw runtime_error("setting has not been terminated");
	}
}

}
/* namespace siena */
