#include "phase2_ce.h"
#include "../utils/utils.h" // For any_int, rng, d_set_to_str
#include <cmath>
#include <limits>
#include <iostream>
#include <map>
#include <vector>
#include <set>

// --- Helper Functions for Locating Algorithm ---
// (These are unchanged)

/**
 * @brief Checks if a given row covers a d-set.
 * A row covers a d-set if it covers all interactions within that d-set.
 */
bool row_covers_d_set(const std::vector<v_type>& row, const d_set_type& d_set) {
    if (d_set.empty()) {
        return false; 
    }
    for (const auto& interaction : d_set) {
        // Check if row covers this single interaction
        bool covers_interaction = true;
        if (interaction.first.empty()) {
             continue;
        }
        for (size_t i = 0; i < interaction.first.size(); ++i) {
            // Bounds check
            if (interaction.first[i] >= row.size() || interaction.second.size() <= i) {
                covers_interaction = false; 
                break;
            }
            if (row[interaction.first[i]] != interaction.second[i]) {
                covers_interaction = false;
                break;
            }
        }
        // If the row fails to cover even one interaction, it fails to cover the d-set
        if (!covers_interaction) {
            return false;
        }
    }
    // If we get here, the row covered all interactions in the d-set
    return true;
}

/**
 * @brief Generates a single random row based on the levels in 'vs'.
 */
std::vector<v_type> generate_random_row(int k, const vs_type& vs) {
    std::vector<v_type> row(k);
    for(int i = 0; i < k; ++i) {
        row[i] = any_int(rng) % vs[i];
    }
    return row;
}


// --- Main Greedy Algorithm Function ---

void run_phase_2_ce(LocatingArray *array) {
	
	int k = array->k;
    const vs_type& vs = array->vs; // Get vs
	int M = array->array.size(); // Start counting from existing rows

    // A set of indices into array->undistinguished_pairs
    std::set<size_t> pair_indices;
    for(size_t i = 0; i < array->undistinguished_pairs.size(); ++i) {
        pair_indices.insert(i);
    }

    if (!pair_indices.empty()) {
        std::cout << "  (CE) DEBUG: Trying to distinguish " << pair_indices.size() << " pairs. Example:" << std::endl;
        const auto& [d_set1, d_set2, needed] = array->undistinguished_pairs[*pair_indices.begin()];
        // Use the utility function from utils/utils.h (which must be linked)
        std::cout << "       Pair 0: " << d_set_to_str(d_set1) << " vs " << d_set_to_str(d_set2) << " (already separated " << needed << " times)" << std::endl;
    }

    // --- NEW: Using greedy, column-by-column construction ---
    // This implements the "Conditional Expectation" heuristic.
    std::cout << "  (CE) Using greedy column-by-column strategy." << std::endl;
    
    // Number of random samples to take when evaluating each value
    const int SAMPLES = 10; 


	while (!pair_indices.empty()) {
		M++;

        // --- Greedy Row Construction (Unchanged) ---
        std::vector<v_type> currentRow(k); // We will build this row greedily
        
        for(int j = 0; j < k; ++j) { // For each column j
            v_type best_v_for_col = 0;     // Best value for this column
            int best_score_for_col = -1; // Best score seen for this column

            // Try every possible value 'v' for column 'j'
            for(v_type v = 0; v < vs[j]; ++v) {
                currentRow[j] = v; // Set the value for this column
                int score_for_v = 0;

                // To score this choice, we randomly fill the *rest* of the
                // row SAMPLES times and sum the scores.
                for (int s = 0; s < SAMPLES; ++s) {
                    
                    // Fill columns j+1 to k-1 randomly
                    for (int r = j + 1; r < k; ++r) {
                        currentRow[r] = any_int(rng) % vs[r];
                    }

                    // Now that `currentRow` is complete, score it
                    // against all remaining pairs
                    for (const auto& index : pair_indices) {
                        const auto& [d_set1, d_set2, needed] = array->undistinguished_pairs[index];
                        
                        bool covers1 = row_covers_d_set(currentRow, d_set1);
                        bool covers2 = row_covers_d_set(currentRow, d_set2);

                        if ((covers1 && !covers2) || (!covers1 && covers2)) {
                            score_for_v++;
                        }
                    }
                } // End sampling loop

                // If this value 'v' gave the best score so far, keep it
                if (score_for_v > best_score_for_col) {
                    best_score_for_col = score_for_v;
                    best_v_for_col = v;
                }
            } // End value loop

            // We've tried all values for col j. Lock in the best one.
            currentRow[j] = best_v_for_col;
        }
        // --- END GREEDY HEURISTIC ---

        // `currentRow` is now the complete, greedily-constructed row
        std::vector<v_type> bestRow = currentRow;

        // Calculate the *actual* score for this row for logging
        int bestScore = 0;
        for (const auto& index : pair_indices) {
            const auto& [d_set1, d_set2, needed] = array->undistinguished_pairs[index];
            bool covers1 = row_covers_d_set(bestRow, d_set1);
            bool covers2 = row_covers_d_set(bestRow, d_set2);
            if ((covers1 && !covers2) || (!covers1 && covers2)) {
                bestScore++;
            }
        }

		// --- ADD THE *BEST* ROW TO THE OBJECT'S ARRAY ---
		array->array.push_back(bestRow);

		// --- *** LOGIC FIX HERE *** ---
        // --- UPDATE THE LIST OF UNDISTINGUISHED PAIRS ---
		std::vector<size_t> distinguished_this_round;

        for (const auto& index : pair_indices) {
            // Get the tuple by reference
            auto& pair_tuple = array->undistinguished_pairs[index];
            const auto& d_set1 = std::get<0>(pair_tuple);
            const auto& d_set2 = std::get<1>(pair_tuple);
            
            // This int is 'times_separated_already', not 'needed'
            int& times_separated = std::get<2>(pair_tuple); 
            
            bool covers1 = row_covers_d_set(bestRow, d_set1);
            bool covers2 = row_covers_d_set(bestRow, d_set2);

            // Check for distinguishing
            if ((covers1 && !covers2) || (!covers1 && covers2)) {
                times_separated++; // <-- INCREMENT the separation count
                
                // Check if it has now met the lambda requirement
                if (times_separated >= array->lambda) { 
                    distinguished_this_round.push_back(index);
                }
            }
        } 
		
        // Remove the pairs that are now fully distinguished
		for (const auto& index : distinguished_this_round) {
			pair_indices.erase(index);
		}
        // --- *** END LOGIC FIX *** ---
		
		if (M % 10 == 0 || pair_indices.empty()) {
			std::cout << "  (CE) Row " << M << " built (best score: " << bestScore << "). "
					  << pair_indices.size() << " pairs remaining to distinguish." << std::endl;
		}

        // --- SAFETY BREAK ---
        if (M > (int(array->array.size()) + k * 20) && M > 200) { 
             std::cout << "  (CE) WARNING: Algorithm seems stuck. Forcefully exiting loop." << std::endl;
             std::cout << "  (CE) " << pair_indices.size() << " pairs were left undistinguished." << std::endl;
             break;
        }
	} 
}