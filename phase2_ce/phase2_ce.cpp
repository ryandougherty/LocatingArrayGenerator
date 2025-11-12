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
	
	int remainingRows = N - M - 1; 

	for (auto const& [interaction, needed] : array->uncovered_interactions) {
		double correct_p_base = 1.0;
		for (auto col_idx : interaction.first) {
			// Check for valid level to prevent division by zero
			if (array->vs[col_idx] > 0) { 
				correct_p_base *= (1.0 / array->vs[col_idx]);
			} else {
				correct_p_base = 0.0; // Cannot be covered if level is 0
				break;
			}
		}
		bool conflicts = false;
        // This must be calculated per-interaction for variable levels
        double probThisRowWillCover = 1.0; 

		for (size_t i = 0; i < interaction.first.size(); ++i) {
			int col = interaction.first[i];
			int val = interaction.second[i];

			if (partialRow[col] != -1) { // If this column is already fixed in the partial row
				if (partialRow[col] != val) {
					conflicts = true; 
                    probThisRowWillCover = 0.0; // This row can't cover it
					break;
				}
                // If partialRow[col] == val, this part matches. probThisRowWillCover stays (1.0 * 1.0)
			} else { 
				// This column is unfixed. The probability of hitting it is 1 / (levels for this col)
                probThisRowWillCover *= (1.0 / array->vs[col]);
			}
		}

		double expectedUncovered;
		if (conflicts) {
            // This row conflicts with the interaction, so it can't cover it.
            // probThisRowWillCover is 0.0.
			expectedUncovered = getUncoverProb(needed, remainingRows, correct_p_base);
		} else {
            // probThisRowWillCover now holds the correct probability (e.g., 1/2 * 1/3 = 1/6)
			expectedUncovered =
				probThisRowWillCover * getUncoverProb(needed - 1, remainingRows, correct_p_base) +
				(1.0 - probThisRowWillCover) * getUncoverProb(needed, remainingRows, correct_p_base);
		}

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
		for (int col = 0; col < k; ++col) {
			double bestExpected = std::numeric_limits<double>::max();
			int bestVal = 0;

			for (int val = 0; val < array->vs[col]; ++val) {
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
			for (size_t i = 0; i < inter.first.size(); ++i) {
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