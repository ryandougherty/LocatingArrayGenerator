/* ----------------------------------------------------------------------------
 * phase2_greedy.cpp
 *
 * MODIFIED: Uses d_set_id from InteractionCodec. Decodes on-the-fly for
 * row coverage checks.
 * ----------------------------------------------------------------------------
 */

#include "phase2_greedy.h"
#include "../utils/utils.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <map>
#include <vector>
#include <set>

// --- Helper Functions ---

/**
 * @brief Checks if a given row covers a d-set.
 * A row covers a d-set if it covers at least one interaction within that d-set.
 * (Because R(D) = UNION of R(I) for each I in D.)
 */
static bool row_covers_d_set(const std::vector<v_type>& row, const d_set_type& d_set) {
    if (d_set.empty()) {
        return false; 
    }
    for (const auto& interaction : d_set) {
        bool covers_interaction = true;
        if (interaction.first.empty()) {
             continue;
        }
        for (size_t i = 0; i < interaction.first.size(); ++i) {
            if (interaction.first[i] >= row.size() || interaction.second.size() <= i) {
                covers_interaction = false; 
                break;
            }
            if (row[interaction.first[i]] != interaction.second[i]) {
                covers_interaction = false;
                break;
            }
        }
        // ANY semantics: if this interaction is covered, the d-set is covered
        if (covers_interaction) {
            return true;
        }
    }
    // No interaction was covered
    return false;
}

/**
 * @brief Generates a single random row based on the levels in 'vs'.
 */
static std::vector<v_type> generate_random_row(int k, const vs_type& vs) {
    std::vector<v_type> row(k);
    for(int i = 0; i < k; ++i) {
        row[i] = any_int(rng) % vs[i];
    }
    return row;
}


// --- Main Greedy Algorithm Function ---

void run_phase_2_greedy(LocatingArray *array) {
	
	int k = array->k;
    const vs_type& vs = array->vs;
	int M = array->array.size();
    const InteractionCodec& codec = array->codec;
    bool is_detecting = array->is_detecting;

    // For detecting: pair is (X_singleton, T_dset). A row separates if it covers X but not T.
    // For locating:  pair is (D1, D2). A row separates if it covers one but not the other.
    auto row_separates = [&](bool covers1, bool covers2) -> bool {
        if (is_detecting)
            return covers1 && !covers2;
        else
            return covers1 != covers2;
    };

    // A set of indices into array->undistinguished_pairs
    std::set<size_t> pair_indices;
    for(size_t i = 0; i < array->undistinguished_pairs.size(); ++i) {
        pair_indices.insert(i);
    }

    if (!pair_indices.empty()) {
        std::cout << "  (CE) DEBUG: Trying to distinguish " << pair_indices.size() << " pairs. Example:" << std::endl;
        const auto& [id1, id2, needed] = array->undistinguished_pairs[*pair_indices.begin()];
        // Decode for display
        d_set_type d_set1 = codec.decode_d_set(id1);
        d_set_type d_set2 = codec.decode_d_set(id2);
        std::cout << "       Pair 0: " << d_set_to_str(d_set1) << " vs " << d_set_to_str(d_set2) << " (already separated " << needed << " times)" << std::endl;
    }

    std::cout << "  (CE) Using greedy column-by-column strategy." << std::endl;
    
    const int SAMPLES = 10; 
    int consecutive_zero_score = 0;
    bool using_random_fallback = false;
    int random_fail_attempts = 0;

    // Cache: d_set_id -> decoded d_set_type (avoids repeated decoding)
    std::unordered_map<d_set_id, d_set_type> decode_cache;
    auto get_decoded = [&](d_set_id id) -> const d_set_type& {
        auto it = decode_cache.find(id);
        if (it != decode_cache.end()) return it->second;
        decode_cache[id] = codec.decode_d_set(id);
        return decode_cache[id];
    };

	while (!pair_indices.empty()) {
		M++;

        std::vector<v_type> bestRow(k);

        if (!using_random_fallback) {
            // --- Greedy column-by-column construction ---
            std::vector<v_type> currentRow(k);
            
            for(int j = 0; j < k; ++j) {
                v_type best_v_for_col = 0;
                int best_score_for_col = -1;

                for(v_type v = 0; v < vs[j]; ++v) {
                    currentRow[j] = v;
                    int score_for_v = 0;

                    for (int s = 0; s < SAMPLES; ++s) {
                        
                        for (int r = j + 1; r < k; ++r) {
                            currentRow[r] = any_int(rng) % vs[r];
                        }

                        for (const auto& index : pair_indices) {
                            const auto& [id1, id2, needed] = array->undistinguished_pairs[index];
                            
                            const d_set_type& d_set1 = get_decoded(id1);
                            const d_set_type& d_set2 = get_decoded(id2);

                            bool covers1 = row_covers_d_set(currentRow, d_set1);
                            bool covers2 = row_covers_d_set(currentRow, d_set2);

                            if (row_separates(covers1, covers2)) {
                                score_for_v++;
                            }
                        }
                    }

                    if (score_for_v > best_score_for_col) {
                        best_score_for_col = score_for_v;
                        best_v_for_col = v;
                    }
                }

                currentRow[j] = best_v_for_col;
            }

            bestRow = currentRow;
        } else {
            // --- Random row fallback ---
            // The greedy has plateaued; just generate a random row.
            // The diagnostic proved random rows always fix remaining pairs.
            for (int c = 0; c < k; ++c) {
                bestRow[c] = any_int(rng) % vs[c];
            }
        }

        // Score the row
        int bestScore = 0;
        for (const auto& index : pair_indices) {
            const auto& [id1, id2, needed] = array->undistinguished_pairs[index];
            const d_set_type& d_set1 = get_decoded(id1);
            const d_set_type& d_set2 = get_decoded(id2);
            bool covers1 = row_covers_d_set(bestRow, d_set1);
            bool covers2 = row_covers_d_set(bestRow, d_set2);
            if (row_separates(covers1, covers2)) {
                bestScore++;
            }
        }

        // Detect plateau and switch to random fallback
        if (bestScore == 0) {
            consecutive_zero_score++;
            if (!using_random_fallback && consecutive_zero_score >= 3) {
                std::cout << "  (CE) Greedy plateaued (3 consecutive zero-score rows). "
                          << "Switching to random row generation.\n";
                using_random_fallback = true;
                // Don't add this useless row, retry with random
                M--;
                continue;
            }
            if (using_random_fallback) {
                // Random row didn't help this time, just skip it
                M--;
                random_fail_attempts++;
                if (random_fail_attempts > 10000) {
                    std::cout << "  (CE) WARNING: 10000 random rows failed. "
                              << pair_indices.size() << " pairs left.\n";
                    break;
                }
                continue;
            }
        } else {
            consecutive_zero_score = 0;
        }

		array->array.push_back(bestRow);

		std::vector<size_t> distinguished_this_round;

        for (const auto& index : pair_indices) {
            auto& pair_tuple = array->undistinguished_pairs[index];
            const d_set_id id1 = std::get<0>(pair_tuple);
            const d_set_id id2 = std::get<1>(pair_tuple);
            int& times_separated = std::get<2>(pair_tuple); 
            
            const d_set_type& d_set1 = get_decoded(id1);
            const d_set_type& d_set2 = get_decoded(id2);

            bool covers1 = row_covers_d_set(bestRow, d_set1);
            bool covers2 = row_covers_d_set(bestRow, d_set2);

            if (row_separates(covers1, covers2)) {
                times_separated++;
                
                if (times_separated >= array->lambda) { 
                    distinguished_this_round.push_back(index);
                }
            }
        } 
		
		for (const auto& index : distinguished_this_round) {
			pair_indices.erase(index);
		}
		
		if (M % 10 == 0 || pair_indices.empty()) {
			std::cout << "  (CE) Row " << M << " built"
                      << (using_random_fallback ? " [random]" : "")
                      << " (best score: " << bestScore << "). "
					  << pair_indices.size() << " pairs remaining to distinguish." << std::endl;
		}

        if (M > (int(array->array.size()) + k * 20) && M > 200) { 
             std::cout << "  (CE) WARNING: Algorithm seems stuck. Forcefully exiting loop." << std::endl;
             std::cout << "  (CE) " << pair_indices.size() << " pairs were left undistinguished." << std::endl;
             break;
        }
	} 
}
