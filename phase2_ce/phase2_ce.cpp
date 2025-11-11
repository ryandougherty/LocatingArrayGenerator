#include "phase2_ce.h"
#include "../utils/utils.h" // For math::combinations
#include <cmath>
#include <limits>
#include <iostream>
#include <map>
#include <vector>

// --- Helper Functions for CE Algorithm ---

double getUncoverProb(int needed, int remainingRows, double p) {
	if (needed <= 0) {
		return 0.0; 
	}
    if (remainingRows <= 0) {
        return 1.0; 
    }

	double probOfFailure = 0.0;
	double q = 1.0 - p;

	for (int i = 0; i < needed; ++i) {
        if (i > remainingRows) break; 
		// --- USE THE NAMESPACED FUNCTION ---
		probOfFailure += math::combinations(remainingRows, i) * std::pow(p, i) * std::pow(q, remainingRows - i);
	}
	return probOfFailure;
}

double calculateTotalExpectedUncovered(LocatingArray *array, const std::vector<v_type>& partialRow, int M, int N) {
	double totalExpected = 0.0;
	
	double p_base = std::pow(1.0 / array->v, array->t);
	int remainingRows = N - M - 1; 

	// --- USE THE OBJECT'S MAP ---
	for (auto const& [interaction, needed] : array->uncovered_interactions) {
		
		int unfixedColsInT = 0;
		bool conflicts = false;

		// --- USE interaction.cols and interaction.vals (from your struct) ---
		for (size_t i = 0; i < interaction.first.size(); ++i) {
			int col = interaction.first[i];
			int val = interaction.second[i];

			if (partialRow[col] != -1) { 
				if (partialRow[col] != val) {
					conflicts = true; 
					break;
				}
			} else { 
				unfixedColsInT++;
			}
		}

		double probThisRowWillCover;
		if (conflicts) {
			probThisRowWillCover = 0.0;
		} else {
			probThisRowWillCover = std::pow(1.0 / array->v, unfixedColsInT);
		}

		double expectedUncovered =
			probThisRowWillCover * getUncoverProb(needed - 1, remainingRows, p_base) +
			(1.0 - probThisRowWillCover) * getUncoverProb(needed, remainingRows, p_base);

		totalExpected += expectedUncovered;
	}
	return totalExpected;
}

int calculateInitialN(int k, int v, int t, int lambda) {
	double N_double = (lambda * std::pow(v, t)) * (1.0 + std::log(k));
    if (N_double < 10) N_double = 10;
    if (N_double > 50000) N_double = 50000; 
	return static_cast<int>(N_double);
}


// --- Main CE Algorithm Function ---

void run_phase_2_ce(LocatingArray *array) {
	
	int k = array->k;
	int v = array->v;
	int t = array->t;
	int lambda = array->lambda;

	int N = calculateInitialN(k, v, t, lambda);
	std::cout << "  (CE) Target N heuristic: " << N << std::endl;

	int M = 0; 

	while (!array->uncovered_interactions.empty()) {
		M++;
		std::vector<v_type> newRow(k, -1); // -1 = "unfixed"
        for (int i = 0; i < k; ++i) {
            // If 'vs' is available on 'array', use that. Otherwise, use 'v'.
            newRow[i] = rand() % (array->vs.empty() ? array->v : array->vs[i]); 
        }
		for (int col = 0; col < k; ++col) {
			double bestExpected = std::numeric_limits<double>::max();
			int bestVal = 0;

			for (int val = 0; val < v; ++val) {
				newRow[col] = val; 
				double expected = calculateTotalExpectedUncovered(array, newRow, M - 1, N);
				if (expected < bestExpected) {
					bestExpected = expected;
					bestVal = val;
				}
			}
			newRow[col] = bestVal;
		}

		// --- ADD THE ROW TO THE OBJECT'S ARRAY ---
		array->array.push_back(newRow);

		// --- UPDATE THE OBJECT'S MAP ---
		std::vector<interaction_type> covered_this_round;
        // Note: We iterate with a reference to the map
		for (auto it = array->uncovered_interactions.begin(); it != array->uncovered_interactions.end(); ++it) {
			// const interaction_type& interaction_str = it->first;
            
            const interaction_type& inter = it->first; // ASSUME key is the 'interaction' struct
			
			bool covers = true;
			for (int i = 0; i < t; ++i) {
				if (newRow[inter.first[i]] != inter.second[i]) {
					covers = false;
					break;
				}
			}

			if (covers) {
				it->second--; 
				if (it->second == 0) {
					covered_this_round.push_back(inter);
				}
			}
		}

		for (const interaction_type& inter : covered_this_round) {
			array->uncovered_interactions.erase(inter);
		}
		
		if (M % 10 == 0 || array->uncovered_interactions.empty()) {
			std::cout << "  (CE) Row " << M << " built. "
					  << array->uncovered_interactions.size() << " interactions remaining to cover." << std::endl;
		}
	}
}