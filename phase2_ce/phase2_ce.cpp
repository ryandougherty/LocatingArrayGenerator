#include "phase2_ce.h"
#include "../utils/utils.h" // For any_int, rng, d_set_to_str
#include <cmath>
#include <limits>
#include <iostream>
#include <map>
#include <vector>
#include <set>

// --- Helper Functions for Locating Algorithm ---

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
        std::cout << "       Pair 0: " << d_set_to_str(d_set1) << " vs " << d_set_to_str(d_set2) << std::endl;
    }

    // --- NEW: Set number of candidates to generate per new row ---
    // Let's try k*2 or 100, whichever is larger.
    int num_candidates = std::max(100, k * 2);
    std::cout << "  (CE) Using randomized greedy strategy with " << num_candidates << " candidates per row." << std::endl;


	while (!pair_indices.empty()) {
		M++;

        // --- NEW: Randomized Greedy Heuristic ---
        std::vector<v_type> bestRow = generate_random_row(k, array->vs);
        int bestScore = -1;

        // Generate and score N candidate rows
        for(int i = 0; i < num_candidates; ++i) {
            std::vector<v_type> currentRow = generate_random_row(k, array->vs);
            int currentScore = 0;

            // Score this row by checking how many *remaining* pairs it distinguishes
            for (const auto& index : pair_indices) {
                const auto& [d_set1, d_set2, needed] = array->undistinguished_pairs[index];
                
                bool covers1 = row_covers_d_set(currentRow, d_set1);
                bool covers2 = row_covers_d_set(currentRow, d_set2);

                if ((covers1 && !covers2) || (!covers1 && covers2)) {
                    currentScore++;
                }
            }

            if (currentScore > bestScore) {
                bestScore = currentScore;
                bestRow = currentRow;
            }
        }
        // --- END NEW HEURISTIC ---

		// --- ADD THE *BEST* ROW TO THE OBJECT'S ARRAY ---
		array->array.push_back(bestRow);

		// --- UPDATE THE LIST OF UNDISTINGUISHED PAIRS ---
		std::vector<size_t> distinguished_this_round;

        for (const auto& index : pair_indices) {
            // Get the tuple by reference
            auto& pair_tuple = array->undistinguished_pairs[index];
            const auto& d_set1 = std::get<0>(pair_tuple);
            const auto& d_set2 = std::get<1>(pair_tuple);
            int& needed = std::get<2>(pair_tuple);
            
            bool covers1 = row_covers_d_set(bestRow, d_set1);
            bool covers2 = row_covers_d_set(bestRow, d_set2);

            // Check for distinguishing
            if ((covers1 && !covers2) || (!covers1 && covers2)) {
                needed--; 
                if (needed == 0) {
                    distinguished_this_round.push_back(index);
                }
            }
        } 
		
        // Remove the pairs that are now fully distinguished
		for (const auto& index : distinguished_this_round) {
			pair_indices.erase(index);
		}
		
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