/*
 * ComposableSetting.cpp
 *
 *  Created on: May 5, 2015
 *      Author: ortmann
 */

#include <stdexcept>
#include <model/settings/ComposableSetting.h>

#include "../../network/iterators/UnionTieIterator.h"
#include "../../network/iterators/SingleIterator.h"

using namespace std;

namespace siena {

ComposableSetting::~ComposableSetting() {
	delete lpfirstSetting;
	delete lpsecondSetting;
}

ComposableSetting::ComposableSetting(Setting* firstSetting,
		Setting* secondSetting) :
		Setting(), //
		lpfirstSetting(firstSetting), //
		lpsecondSetting(secondSetting), //
		lpSteps(0), //
		lpPermittedSteps(0) {
}

void ComposableSetting::initSetting(Network* const lpNetwork) {
	lpfirstSetting->initSetting(lpNetwork);
	lpsecondSetting->initSetting(lpNetwork);
}

void ComposableSetting::terminateSetting(Network* const lpNetwork) {
	lpfirstSetting->terminateSetting(lpNetwork);
	lpsecondSetting->terminateSetting(lpNetwork);
}

void ComposableSetting::initDyadicSetting(const std::map<int, double>& row,
		int ego) {
	lpfirstSetting->initDyadicSetting(row, ego);
	lpsecondSetting->initDyadicSetting(row, ego);
}

void ComposableSetting::initPermittedSteps(const bool* const permitted) {
	if (lpPermittedSteps == 0) {
	lpfirstSetting->initPermittedSteps(permitted);
	lpsecondSetting->initPermittedSteps(permitted);
	ITieIterator* iter1 = lpfirstSetting->getPermittedSteps();
	ITieIterator* iter2 = lpsecondSetting->getPermittedSteps();
		SingleIterator sIter(ego());
		UnionTieIterator uIter(*iter1, sIter);
		lpPermittedSteps = new UnionTieIterator(uIter, *iter2);
	delete iter1;
	delete iter2;
	} else {		
	throw runtime_error("setting has not been terminated");
	}
}

ITieIterator* ComposableSetting::getSteps() {
	if (lpSteps != 0) {
		return lpSteps->clone();
	}
	throw runtime_error("setting has not been initialized");
}

int ComposableSetting::getSize() {
	return lpSteps->size();
}

ITieIterator* ComposableSetting::getPermittedSteps() {
	return lpPermittedSteps->clone();
}

int ComposableSetting::getPermittedSize() {
	return lpPermittedSteps->size();
}

void ComposableSetting::initSetting() {
	lpfirstSetting->initSetting(ego());
	lpsecondSetting->initSetting(ego());

	ITieIterator* iter1 = lpfirstSetting->getSteps();
	ITieIterator* iter2 = lpsecondSetting->getSteps();

	SingleIterator sIter(ego());
	UnionTieIterator uIter(*iter1, sIter);

	lpSteps = new UnionTieIterator(uIter, *iter2);
	delete iter1;
	delete iter2;
}

void ComposableSetting::terminateSetting() {
	if (lpSteps != 0) {
	lpfirstSetting->terminateSetting(ego());
	lpsecondSetting->terminateSetting(ego());
	delete lpPermittedSteps;
	delete lpSteps;
	lpPermittedSteps = 0;
	lpSteps = 0;
	} else {
		throw runtime_error("setting has not been initialized");
	}
}

}
/* namespace siena */
