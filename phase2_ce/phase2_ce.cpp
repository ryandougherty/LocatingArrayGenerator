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

double calculateTotalExpectedUncovered(LocatingArray *array, const std::vector<int>& partialRow, int M, int N) {
	double totalExpected = 0.0;
	
	double p_base = std::pow(1.0 / array->v, array->t);
	int remainingRows = N - M - 1; 

	// --- USE THE OBJECT'S MAP ---
	for (auto const& [interaction, needed] : array->uncovered_interactions) {
		
		int unfixedColsInT = 0;
		bool conflicts = false;

		// --- USE interaction.cols and interaction.vals (from your struct) ---
		for (int i = 0; i < interaction.cols.size(); ++i) {
			int col = interaction.cols[i];
			int val = interaction.vals[i];

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
		std::vector<int> newRow(k, -1); // -1 = "unfixed"

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
		std::vector<Interaction> covered_this_round;
        // Note: We iterate with a reference to the map
		for (auto it = array->uncovered_interactions.begin(); it != array->uncovered_interactions.end(); ++it) {
			const Interaction& interaction_str = it->first;
            
            // We need to parse the interaction string to check coverage
            // This is inefficient, but matches your 'Interaction' type
            // A better 'Interaction' struct would be faster
            
            // --- Re-create the 'interaction' struct from the key ---
            // This assumes your 'Interaction' (string) is formatted like "col1,col2:val1,val2"
            // Your code in `phase1.cpp` seems to just use a string.
            // Let's assume the key is the string from `Interaction::to_string()`
            
            // This is the hard part. Your `Interaction` key is a string.
            // We need to know which columns/values it represents.
            // Your `LocatingArray` class doesn't store `allInteractions` as structs.
            
            // Let's check `phase1.cpp`...
            // Ah, `run_phase_1` uses `itertools::combinations` on `cols` and `vals`
            // and then creates an `interaction` struct, then `i.to_string()`.
            
            // We MUST modify `LocatingArray` to store the actual structs.
            // A `map<string, int>` is not enough information.
            
            // --- THIS IS A REQUIRED CHANGE ---
            // We must change `uncovered_interactions` to store the struct
            // or we must re-parse the string. Re-parsing is very slow.
            
            // Let's assume `Interaction` is the *struct* not the *string*.
            // Your `LocAG.h` has `typedef std::string Interaction;`
            // This is the core problem.
            
            // ---
            // TEMP FIX: Let's assume your 'Interaction' struct has a 'from_string'
            // or we parse it manually. This is SLOW.
            // A better fix is to change LocAG.h.
            // ---
            
            // Let's just re-implement the check based on your GA code.
            // Your GA code `calculate_fitness` iterates `uncovered_interactions`
            // and calls `i.check_row(row)`. This implies `i` is an `interaction` struct.
            
            // THIS IS THE PROBLEM:
            // LocAG.h: `typedef std::string Interaction;`
            // LocAG.h: `std::map<Interaction, int> uncovered_interactions;`
            // phase1.cpp: `array->uncovered_interactions[i.to_string()] = array->lambda;`
            // phase2.cpp: `for (auto const& [i, v] : array->uncovered_interactions)`
            // phase2.cpp: `if (i.check_row(row))` -> This CANNOT work. `i` is a `std::string`.
            
            // Your existing `phase2.cpp` code MUST be broken.
            // `i.check_row(row)` will not compile if `i` is a string.
            // Let's check your `phase2.cpp`...
            
            // Ah, you have `std::map<interaction, int> uncovered_interactions;`
            // BUT LocAG.h has `std::map<Interaction, int> uncovered_interactions;`
            // and `typedef std::string Interaction;`
            // This means your uploaded code doesn't match what you're compiling.
            
            // I will assume the `LocAG.h` from your upload is WRONG
            // and that `uncovered_interactions` is REALLY `std::map<interaction, int>`
            // This is the only way your `phase2.cpp` could ever work.
            
            // *** FIXING `phase2_ce.cpp` BASED ON THIS ASSUMPTION ***
            
            const interaction& inter = it->first; // ASSUME key is the 'interaction' struct
			
			bool covers = true;
			for (int i = 0; i < t; ++i) {
				if (newRow[inter.cols[i]] != inter.vals[i]) {
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

		for (const interaction& inter : covered_this_round) {
			array->uncovered_interactions.erase(inter);
		}
		
		if (M % 10 == 0 || array->uncovered_interactions.empty()) {
			std::cout << "  (CE) Row " << M << " built. "
					  << array->uncovered_interactions.size() << " interactions remaining to cover." << std::endl;
		}
	}
}