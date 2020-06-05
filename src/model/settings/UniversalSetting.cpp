/*
 * UniversalSetting.cpp
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#include <stdexcept>
#include "UniversalSetting.h"
#include "../../network/Network.h"
#include "../../network/iterators/IntVecIterator.h"


using namespace std;

namespace siena {

UniversalSetting::UniversalSetting() :
		GeneralSetting(), //
		rSteps() {
}

UniversalSetting::~UniversalSetting() {
}

void UniversalSetting::initSetting(Network* const lpNetwork) {
	if (!rSteps.empty()) {
		throw runtime_error("setting has not been terminated");
	}
	rSteps.reserve(lpNetwork->m());
	for (int i = 0; i < lpNetwork->m(); i++) {
		rSteps.push_back(i);
	}
}

void UniversalSetting::terminateSetting(Network* const lpNetwork) {
	rSteps.clear();
}

void UniversalSetting::initSetting() {
	return;
}

ITieIterator* UniversalSetting::getSteps() {
	if (!rSteps.empty()) {
		return new IntVecIterator(rSteps.begin(), rSteps.end());
	}
	throw runtime_error("setting has not been initialized");
}

int UniversalSetting::getSize() {
	return rSteps.size();
}

} /* namespace siena */
