/*
 * PrimarySetting.cpp
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#undef length
#include <Rinternals.h>
#include "PrimarySetting.h"

#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/iterators/ITieIterator.h"
#include "network/iterators/UnionTieIterator.h"
#include "network/iterators/FilteredIterator.h"
#include "network/iterators/SingleIterator.h"
#include "network/IncidentTieIterator.h"

// #include <algorithm> // std::find
// #include "siena07utilities.h"


namespace siena {

PrimarySetting::PrimarySetting() :
		GeneralSetting(),
		lpNetwork(0),
		// rDistTwoLayer(),
		// lpiter(0),
		lPLayer() {
}

PrimarySetting::~PrimarySetting() {
	// if (lpiter != 0) {
	// 	delete lpiter;
	// }
}

// static SEXP asInt(int i) {
// 	SEXP ri = allocVector(INTSXP, 1);
// 	INTEGER(ri)[0] = i;
// 	return ri;
// }

// static std::vector<int> str(SEXP exp, SEXP l, int ii) {
// 	SEXP call = PROTECT(allocList(4));
// 	SET_TYPEOF(call, LANGSXP);

// 	SETCAR(call, install("debugg"));
// 	SETCAR(CDR(call), exp);
// 	SETCAR(CDDR(call), l);
// 	SETCAR(CDDDR(call), asInt(ii));

// 	SEXP r = eval(call, R_GlobalEnv);

// 	std::vector<int> ret;
// 	for (int i = 0; i < Rf_length(r); i++)
// 		ret.push_back(INTEGER(r)[i]);

// 	UNPROTECT(1);
// 	return ret;
// }

// static void check(const Network * pLayer, const Network * pCounts, const Network * pNetwork, int ego) {
// 	// all r neigbors
// 	std::vector<int> rmore = str(net_to_sexp(pNetwork), net_to_sexp(pCounts), ego);
// 	int rs = rmore.size();
// 	std::vector<int> rmiss;
// 	// delete items in rmore
// 	ITieIterator* n = pLayer->outTies(ego).clone();
// 	while (n->valid()) {
// 		// find in r list
// 		std::vector<int>::iterator p = std::find(rmore.begin(), rmore.end(), n->actor());
// 		if (p != rmore.end()) {
// 			// if found, delete
// 			rmore.erase(p);
// 		} else {
// 			// if not in r, remeber it
// 			rmiss.push_back(n->actor());
// 		}
// 		n->next();
// 	}
// 	delete n;

// 	if (rmiss.size() > 0 || rmore.size() > 0) {
// 		std::stringstream ss;
// 		ss << "fail: r=" << rs << ",c=" << pLayer->outDegree(ego);
// 		ss << ",not in r=";
// 		for(size_t i = 0; i < rmiss.size(); ++i) ss << rmiss[i] << " ";
// 		ss << ",additional=";
// 		for(size_t i = 0; i < rmore.size(); ++i) ss << rmore[i] << " ";
// 		std::string s = ss.str();
// 		LOGS(Priority::ERROR) << ss.str();
// 		throw "fail";
// 	}
// }

void PrimarySetting::initSetting(Network* const lpNetwork) {
	// lpNetwork->addNetworkChangeListener(&rDistTwoLayer);
	lpNetwork->addNetworkChangeListener(&lPLayer);
	this->lpNetwork = lpNetwork;

	// for (int i = 0; i < lpNetwork->n(); i++) {
	// 	LOGS(Priority::ERROR) << "init setting " << i;
	// 	check(lPLayer.pLayer(), lPLayer.pCounts(), lpNetwork, i);
	// }
}

void PrimarySetting::terminateSetting(Network* const lpNetwork) {
	// lpNetwork->removeNetworkChangeListener(&rDistTwoLayer);
	// rDistTwoLayer.clear(lpNetwork->n());
	lpNetwork->removeNetworkChangeListener(&lPLayer);
	this->lpNetwork = 0;
}

void PrimarySetting::initSetting() {
	// if (lpiter == 0) {
	// 	IncidentTieIterator iter1 = lpNetwork->inTies(ego());
	// 	IncidentTieIterator iter2 = lpNetwork->outTies(ego());
	// 	UnionTieIterator uIter1(iter1, iter2);
	// 	iter1 = rDistTwoLayer.getDistanceTwoNeighbors(ego());
	// 	SingleIterator egoIter(ego());
	// 	UnionTieIterator uIter2(iter1, egoIter);
	// 	lpiter = new UnionTieIterator(uIter1, uIter2);
	// } else {
	// 	LOGS(Priority::ERROR)<< "setting has not been terminated\n";
	// 	throw "setting has not been terminated";
	// }
}

void PrimarySetting::terminateSetting() {
	GeneralSetting::terminateSetting();
	// if (lpiter != 0) {
	// 	delete lpiter;
	// 	lpiter=0;
	// 	GeneralSetting::terminateSetting();
	// } else {
	// 	LOGS(Priority::ERROR)<< "setting has not been initialized\n";
	// 	throw "setting has not been initialized";
	// }
}

ITieIterator* PrimarySetting::getSteps() {
	// if (lpiter == 0) throw "setting has not been initialized";
	// for (int i = 0; i < lpNetwork->n(); i++) {
	// 	// LOGS(Priority::ERROR) << "get Setps check " << i;
	// 	check(lPLayer.pLayer(), lPLayer.pCounts(), lpNetwork, i);
	// }
	return this->pPrimaryNetwork()->outTies(ego()).clone();
	}

int PrimarySetting::getSize() {
	return this->pPrimaryNetwork()->outDegree(ego());
}

const Network * PrimarySetting::pPrimaryNetwork() const {
	return lPLayer.pLayer();
}

} // namespace siena
