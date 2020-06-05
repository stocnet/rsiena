/*
 * MeetingSetting.cpp
 *
 *  Created on: 23.06.2014
 *      Author: ortmann
 */

#include <model/settings/MeetingSetting.h>
#include "../../network/iterators/UnionTieIterator.h"
#include "../../network/iterators/SingleIterator.h"
#include "../../utils/Random.h"
#include "../../network/Network.h"

using namespace std;


namespace siena {

MeetingSetting::~MeetingSetting() {
	delete lpSetting;
	if (lpPermittedSteps != 0) {
		delete lpPermittedSteps;
	}
}

MeetingSetting::MeetingSetting(Setting* setting, const Permission_Type type) :
		Setting(), //
		lpSetting(setting), //
		lpPermittedSteps(0), //
		lPermissionType(type) {
}

void MeetingSetting::initSetting(Network* const lpNetwork) {
	lpSetting->initSetting(lpNetwork);
}

void MeetingSetting::terminateSetting(Network* const lpNetwork) {
	lpSetting->terminateSetting(lpNetwork);
}

void MeetingSetting::initDyadicSetting(const std::map<int, double>& row,
		int ego) {
	lpSetting->initDyadicSetting(row, ego);
}

void MeetingSetting::initPermittedSteps(const bool* const permitted) {
	// THIS WORKS CAUSE WE KNOW THAT PERMITTED STEPS CONTAINS EGO!
	if (lpPermittedSteps == 0) {
	lpSetting->initPermittedSteps(permitted);
		if (lpSetting->getPermittedSize() > 1) {
		ITieIterator* iter = lpSetting->getPermittedSteps();
		if(iter->actor() == ego()) {
			iter->next();
		}
			int pos = nextInt(lpSetting->getPermittedSize() - 1);
		while (pos != 0) {
			iter->next();
			if (iter->actor() != ego()) {
				--pos;
			}
		}
			SingleIterator iter1(ego());
			SingleIterator iter2(iter->actor());
			lpPermittedSteps = new UnionTieIterator(iter1, iter2);
		delete iter;
	} else {
		lpPermittedSteps = new SingleIterator(ego());
	}
	} else {
		throw runtime_error("setting has not been terminated");
	}
}

ITieIterator* MeetingSetting::getSteps() {
	return lpSetting->getSteps();
}

int MeetingSetting::getSize() {
	return lpSetting->getSize();
}

ITieIterator* MeetingSetting::getPermittedSteps() {
	return lpPermittedSteps->clone();
}

int MeetingSetting::getPermittedSize() {
	return lpPermittedSteps->size();
}

void MeetingSetting::initSetting() {
	lpSetting->initSetting(ego());
}

bool MeetingSetting::validate(const Network* const lpNetwork) {
	if (Setting::validate(lpNetwork)) {
	if (lPermissionType != Permission_Type::BOTH) {
			while (lpPermittedSteps->valid()
					&& lpPermittedSteps->actor() == ego()) {
			lpPermittedSteps->next();
		}
		// lpPermittedSteps is since there are at least 2 entries
		// and they are distinct
		bool hasEdge = lpNetwork->hasEdge(ego(), lpPermittedSteps->actor());
		lpPermittedSteps->reset();
		// if up and has edge return false
		if (lPermissionType == Permission_Type::UP && hasEdge) {
			return false;
		}
		// if down and edge no existent return false
		if (lPermissionType == Permission_Type::DOWN && !hasEdge) {
			return false;
		}
	}
	return true;
}
	return false;
}

void MeetingSetting::terminateSetting() {
	lpSetting->terminateSetting(ego());
	if (lpPermittedSteps != 0) 
	{
		delete lpPermittedSteps;
		lpPermittedSteps = 0;
	} 
	else 
	{
		throw runtime_error("setting has not been initialized");
	}
}

}
/* namespace siena */
