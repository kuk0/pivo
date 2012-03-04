#ifndef _CAPRARA_HPP_
#define _CAPRARA_HPP_

	// caprara, cip, ecip defines

	// bounds
#define UB 1
#define LB -1

#define ONE_MEDIAN 0
#define ALL_MEDIAN 1

#define LAZY 1
#define NOTLAZY 0

#include <vector>

using namespace std;



/**
 * wrapper for capraras median solver in grappa
 * @param[in] genomes the genomes
 * @param[in] lower_bound lower bound of the median score,
 *	only solutions with a score >= lower_bound are returned,
 *	if the value is 0 or INT_MAX the best solutions are returned
 * @param[in] upper_bound upper_bound of the median score,
 *	only solutions with a score <= upper_bound are returned,
 *	if the value is 0 or INT_MAX the best solutions are returned
 * @param[in] allmed get all medians
 * @param[in] lazy if in the one-median-version the trivial median score is equal to the lower bound
 *	simply take this one (this is the standard behaviour implemented in GRAPPA .. uhh ohh)
 * @param[out] mediancnt number of medians
 * @param[out] medians the medians
 */

void caprara(const vector<vector<int> > &genomes,
	int lower_bound, int upper_bound,
	int allmed, int lazy,
	int &mediancnt, vector<vector<int> > &medians);

#endif
