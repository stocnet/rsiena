/*
 * Setting.cpp
 *
 *  Created on: 23.07.2014
 *      Author: ortmann
 */

#include "Setting.h"

#include "../../network/iterators/UnionTieIterator.h"
#include "../../network/iterators/SingleIterator.h"
#include "../../network/Network.h"

namespace siena {

void Setting::initSetting(int ego) {
	lego = ego;
	initSetting();
}

void Setting::terminateSetting(int ego) {
	terminateSetting();
	lego = -1;
}

bool Setting::validate(const Network* const lpNetwork) {
	return getPermittedSize() > 1;
}

} /* namespace siena */

