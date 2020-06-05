/*
 * DyadicSetting.cpp
 *
 *  Created on: 24.07.2014
 *      Author: ortmann
 */

#include <stdexcept>
#include "DyadicSetting.h"

#include "../../network/iterators/UnionTieIterator.h"
#include "../../network/iterators/SingleIterator.h"
#include "../../network/iterators/IntDoubleMapIterator.h"


using namespace std;

namespace siena {

DyadicSetting::DyadicSetting() :
		GeneralSetting(), //
		lpiter(0) {
}

DyadicSetting::~DyadicSetting() {
	if (lpiter != 0) {
		delete lpiter;
	}
}

void DyadicSetting::initSetting(Network* const lpNetwork) {
	return;
}

void DyadicSetting::terminateSetting(Network* const lpNetwork) {
	return;
}

void DyadicSetting::initSetting() {
	return;
}

int DyadicSetting::getSize() {
	return lpiter->size();
}

ITieIterator* DyadicSetting::getSteps() {
	return lpiter->clone();
}

void DyadicSetting::initDyadicSetting(const std::map<int, double>& row,
		int ego) {
	if (lpiter == 0) {
	if (row.find(ego) == row.end()) {
			IntDoubleMapIterator iter1(row.begin(), row.end());
			SingleIterator iter2(ego);
			lpiter = new UnionTieIterator(iter1, iter2);
	} else {
		lpiter = new IntDoubleMapIterator(row.begin(), row.end());
	}
	} else {
	throw runtime_error("setting has not been terminated or is used in different contexts");
	}
}

void DyadicSetting::terminateSetting() {
	if (lpiter != 0) 
	{
		delete lpiter;
		lpiter = 0;
		GeneralSetting::terminateSetting();
	} 
	else 
	{
		throw runtime_error("setting has not been initialized");
	}
}

}
/* namespace siena */
