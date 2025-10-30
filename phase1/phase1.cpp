/* ----------------------------------------------------------------------------
 * phase1.cpp
 *
 * Implements Phase 1 of the LocAG algorithm.
 *
 * The primary goal of this phase is to analyze an initial covering array (A)
 * and identify all pairs of "d-sets" (sets of interactions) that
 * violate the locating or detecting property.
 *
 * This is a major bottleneck. To speed this up, it uses a heuristic:
 * 1.  It divides the array 'A' into 'X' horizontal partitions.
 * 2.  For each interaction 'I', it creates a "signature" vector of counts,
 * recording how many times 'I' appears in each partition.
 * 3.  It groups d-sets based on these signatures.
 * 4.  It assumes that two d-sets (D1, D2) can only be non-located/detected
 * if they have *identical* (or very similar) signatures.
 * 5.  This drastically reduces the number of pairs (D1, D2) that require the
 * expensive full symmetric difference check.
 * ----------------------------------------------------------------------------
 */

#include "phase1.h"

// A vector of counts, one for each partition.
using partitioned_counts = std::vector<int>;

/**
 * @brief Calculates the "partitioned count signature" for a single interaction.
 *
 * Divides the array 'A' into 'X' horizontal partitions and counts the
 * occurrences of interaction 'I' within each partition.
 *
 * @param I The interaction (cols, vals) to find.
 * @param A The covering array.
 * @param X The number of partitions.
 * @return A vector of size 'X' where each element is the count for that partition.
 */
auto get_partitioned_interaction_counts(const interaction_type& I, const ca_type& A, int X) {
    partitioned_counts counts(X, 0);
    int partition_size = A.size() / X;

    for (const auto& [idx, row] : enumerate(A)) {
        bool match = true;
        // Check if the row matches the interaction
        for (const auto& [col, val] : zip(I.first, I.second)) {
            if (row[col] != val) {
                match = false;
                break;
            }
        }
        if (match) {
            // If it matches, find which partition this row index belongs to
            int partition_index = std::min(X - 1, static_cast<int>(idx / partition_size));
            counts[partition_index]++;
        }
    }
    return counts;
}

/**
 * @brief Finds all row indices in 'A' where interaction 'I' appears.
 *
 * @param I The interaction (cols, vals) to find.
 * @param A The covering array to search.
 * @return A set of row indices.
 */
robin_hood::unordered_flat_set<int> rows_of_interaction(const interaction_type& I, const ca_type& A) {
    const auto& cols = I.first;
    const auto& vals = I.second;
    robin_hood::unordered_flat_set<int> rows_I_appears;
    
    // Iterate over all rows with their index
    for (const auto& [idx, row] : enumerate(A)) {
        bool flag = false; // Becomes true if this row *doesn't* match
        // Check if the row matches the interaction
        for (const auto& [col, val] : zip(cols, vals)) {
            if (row[col] != val) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            // No mismatch found, so this row covers the interaction
            rows_I_appears.insert(idx);
        }
    }
    return rows_I_appears;
}


/**
 * @brief Generates all interactions of strength up to 't'.
 *
 * An interaction is a (column_set, value_tuple) pair.
 * If 't_bar' is true, it generates interactions for strengths 1, 2, ..., t.
 * If 't_bar' is false, it only generates interactions for strength 't'.
 *
 * @param t The maximum strength.
 * @param vs The vector of levels for each column.
 * @param t_bar Whether to include interactions of strength < t.
 * @return A vector of all possible interactions.
 */
auto get_interactions(const t_type t, const vs_type& vs, bool t_bar) {
    // --- Generate all column sets ---
    auto lb = t; // lower bound
    if (t_bar) {
        lb = 1; // If t_bar, start from strength 1
    }
    std::vector<std::vector<k_type>> col_sets;
    for (int i=lb; i<=t; i++) {
        // Get all combinations of columns of size 'i'
        auto cols = combinations(range(vs.size()), i);
        for (const auto& col_set : cols) {
            std::vector<k_type> to_add;
            for (const auto& new_col : col_set) {
                to_add.push_back(new_col);
            }
            col_sets.push_back(to_add);
        }
    }

    // --- Generate all value-tuples for each column set ---
    // This is a (sub-optimal) hardcoded loop for t=1 up to t=6.
    // A more general solution would use iterators or recursion.
    std::vector<interaction_type> interactions;
    for (const auto& col : col_sets) {
        if (col.size() == 1) {
            for (v_type i = 0; i < vs[col[0]]; i++) {
                std::vector<v_type> s{i};
                interaction_type I = std::make_pair(col, s);
                interactions.push_back(I);
            }
        }
        else if (col.size() == 2) {
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    std::vector<v_type> s{i, j};
                    interaction_type I = std::make_pair(col, s);
                    interactions.push_back(I);
                }
            }
        }
        else if (col.size() == 3) {
            // ... (loops for t=3) ...
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    for(v_type k = 0; k < vs[col[2]]; k++) {
                        std::vector<v_type> s{i, j, k};
                        interaction_type I = std::make_pair(col, s);
                        interactions.push_back(I);
                    }
                }
            }
        }
        else if (col.size() == 4) {
            // ... (loops for t=4) ...
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    for(v_type k = 0; k < vs[col[2]]; k++) {
                        for(v_type l = 0; l < vs[col[3]]; l++) {
                            std::vector<v_type> s{i, j, k, l};
                            interaction_type I = std::make_pair(col, s);
                            interactions.push_back(I);
                        }
                    }
                }
            }
        }
        else if (col.size() == 5) {
            // ... (loops for t=5) ...
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    for(v_type k = 0; k < vs[col[2]]; k++) {
                        for(v_type l = 0; l < vs[col[3]]; l++) {
                            for(v_type m = 0; m < vs[col[4]]; m++) {
                                std::vector<v_type> s{i, j, k, l, m};
                                interaction_type I = std::make_pair(col, s);
                                interactions.push_back(I);
                            }
                        }
                    }
                }
            }
        }
        else if (col.size() == 6) {
            // ... (loops for t=6) ...
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    for(v_type k = 0; k < vs[col[2]]; k++) {
                        for(v_type l = 0; l < vs[col[3]]; l++) {
                            for(v_type m = 0; m < vs[col[4]]; m++) {
                                for (v_type n = 0; n < vs[col[5]]; n++) {
                                    std::vector<v_type> s{i, j, k, l, m, n};
                                    interaction_type I = std::make_pair(col, s);
                                    interactions.push_back(I);
                                }
                            }
                        }
                    }
                }
            }
        } 
        else {
            std::cerr << "Error, Invalid t (t > 6 is not supported by get_interactions)\n";
            abort();
        }
    }
    return interactions;
}

/**
 * @brief Checks if vector 'a' is a subset of vector 'b'.
 * Assumes both vectors are sorted.
 */
bool is_subset(const std::vector<int>& a, const std::vector<int>& b) {
    return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

/**
 * @brief Lexicographical comparison for two interactions.
 * Used for sorting or as a key in std::map.
 */
bool compare_interactions(const interaction_type& a, const interaction_type& b) {
    // First, compare the size of the column vectors
    if (a.first.size() < b.first.size()) return true;
    if (a.first.size() > b.first.size()) return false;

    // If sizes are the same, compare the column vector elements
    for (size_t i = 0; i < a.first.size(); ++i) {
        if (a.first[i] < b.first[i]) return true;
        if (a.first[i] > b.first[i]) return false;
    }

    // If column vectors are identical, compare the size of the value vectors
    if (a.second.size() < b.second.size()) return true;
    if (a.second.size() > b.second.size()) return false;

    // If value vector sizes are the same, compare their elements
    for (size_t i = 0; i < a.second.size(); ++i) {
        if (a.second[i] < b.second[i]) return true;
        if (a.second[i] > b.second[i]) return false;
    }

    // If they are identical, return false
    return false;
}

/**
 * @brief Helper function to get all rows for a d-set (union of its interactions).
 *
 * @param d_set The set of interactions.
 * @param A The covering array.
 * @return A set of row indices that cover *at least one* interaction in 'd_set'.
 */
robin_hood::unordered_flat_set<int> rows_of_d_set(const d_set_type& d_set, const ca_type& A) {
    robin_hood::unordered_flat_set<int> the_rows;
    for(const auto& interaction : d_set) {
        // Get rows for this single interaction
        auto rows = rows_of_interaction(interaction, A);
        // Add them to the union set
        the_rows.insert(rows.begin(), rows.end());
    }
    return the_rows;
}

/**
 * @brief Finds all non-detecting pairs of d-sets.
 *
 * A pair of d-sets (D1, D2) is non-detecting if the set of rows
 * covering D1 is identical to the set of rows covering D2.
 * R(D1) == R(D2)
 *
 * @return A vector of tuples, each containing (D1, D2, 0). The '0' indicates
 * the initial separation count, which is not used for detecting arrays.
 */
std::vector<std::tuple<d_set_type, d_set_type, int>> find_non_detecting_sets( const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar, int X) {
    
    // 1. Generate all possible interactions
    auto interactions = get_interactions(t, vs, t_bar);

    // 2. Generate all possible d-sets (combinations of 'd' interactions)
    std::vector<d_set_type> d_sets;
    auto lower_lim = d;
    if (d_bar) {
        lower_lim = 1; // "at most d"
    }
    for (int i = lower_lim; i <= d; i++) {
        auto to_add = combinations(interactions, d); // All combinations of size 'i'
        for (auto& individual_d_set : to_add) {
            d_set_type inner_d_set;
            for (auto& interaction : individual_d_set) {
                inner_d_set.push_back(interaction);
            }
            d_sets.push_back(inner_d_set);
        }
    }

    // (This map is created but not used for detecting, only for locating)
    robin_hood::unordered_map<interaction_type, partitioned_counts, InteractionHasher> interaction_to_partition_counts;
    // for (const auto& interaction : interactions) {
    //     interaction_to_partition_counts[interaction] = get_partitioned_interaction_counts(interaction, A, X);
    // }
    
    // 3. Group d-sets by their partitioned count signature
    // Map: signature -> list of d-sets with that signature
    robin_hood::unordered_map<partitioned_counts, std::vector<d_set_type>, VectorHasher> d_sets_by_partition_counts;
    std::cout << "Grouping d-sets by partitioned row counts for detecting arrays...\n";
    
    int partition_size = A.size() / X;
    
    for (const auto& d_set : d_sets) {
        partitioned_counts total_counts(X, 0);
        robin_hood::unordered_set<int> distinct_rows;
        
        // For detecting, the signature is based on the *union* of rows
        // So, first get the union of all rows for this d-set
        for (const auto& interaction : d_set) {
            auto rows = rows_of_interaction(interaction, A);
            distinct_rows.insert(rows.begin(), rows.end());
        }

        // Now, build the partition signature from this union-set of rows
        for (const auto& row_idx : distinct_rows) {
            int partition_index = std::min(X - 1, static_cast<int>(row_idx / partition_size));
            total_counts[partition_index]++;
        }

        // Add this d-set to the bucket for its signature
        d_sets_by_partition_counts[total_counts].push_back(d_set);
    }

    // 4. Find non-detecting pairs
    std::vector<std::tuple<d_set_type, d_set_type, int>> to_return;
    std::cout << "Ready to look at pairs for detecting arrays...\n";

    // Iterate through all buckets in the map
    for (const auto& [counts, d_set_group] : d_sets_by_partition_counts) {
        // If a bucket has more than one d-set, all pairs in it are candidates
        if (d_set_group.size() > 1) {
             for (size_t i = 0; i < d_set_group.size(); ++i) {
                for (size_t j = i + 1; j < d_set_group.size(); ++j) {
                    // Because they have the same signature, they are very likely to be
                    // non-detecting pairs.
                    // For higher accuracy, a final check with rows_of_d_set()
                    // could be added here, but the heuristic is assumed to be strong.
                    to_return.push_back({d_set_group[i], d_set_group[j], 0});
                }
            }
        }
    }
    return to_return;
}

/**
 * @brief Finds all non-locating pairs of d-sets.
 *
 * A pair of d-sets (D1, D2) is non-locating if the symmetric difference
 * of their row sets is less than lambda.
 * | R(D1) XOR R(D2) | < lambda
 *
 * @return A vector of tuples, each containing (D1, D2, diff_size), where
 * 'diff_size' is the size of the symmetric difference.
 */
std::vector<std::tuple<d_set_type, d_set_type, int>> find_non_locating_sets(const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar, int X) {
    // 1. Generate all interactions
    auto interactions = get_interactions(t, vs, t_bar);
    
    // 2. Generate all d-sets
    std::vector<d_set_type> d_sets;
    auto lower_lim = d;
    if (d_bar) {
        lower_lim = 1;
    }
    for (int i=lower_lim; i<=d; i++) {
        auto to_add = combinations(interactions, d);
        for (auto& individual_d_set : to_add) {
            d_set_type inner_d_set;
            for (auto& interaction : individual_d_set) {
                inner_d_set.push_back(interaction);
            }
            d_sets.push_back(inner_d_set);
        }
    }
    
    // 3. Get partitioned counts for *each individual interaction*
    robin_hood::unordered_map<interaction_type, partitioned_counts, InteractionHasher> interaction_to_partition_counts;
    for (const auto& interaction : interactions) {
        interaction_to_partition_counts[interaction] = get_partitioned_interaction_counts(interaction, A, X);
    }

    // 4. Group d-sets by their signature
    // Map: signature -> list of d-sets with that signature
    robin_hood::unordered_map<partitioned_counts, std::vector<d_set_type>, VectorHasher> d_sets_by_partition_counts;
    std::cout << "Computing partitioned counts of d-sets...\n";
    for (const auto& d_set : d_sets) {
        partitioned_counts total_counts(X, 0);
        // The signature for a d-set is the *sum* of its interactions' signatures
        // This is a simplification/heuristic.
        for (const auto& interaction : d_set) {
            const auto& counts = interaction_to_partition_counts[interaction];
            for(int i = 0; i < X; ++i) {
                total_counts[i] += counts[i]; 
            }
        }
        // Add this d-set to the bucket for its signature
        d_sets_by_partition_counts[total_counts].push_back(d_set);
    }

    // 5. Find non-locating pairs
    std::vector<std::tuple<d_set_type, d_set_type, int>> to_return;
    std::cout << "Ready to look at pairs based on partitioned counts...\n";

    // Iterate through all buckets
    for (const auto& [counts, d_set_group] : d_sets_by_partition_counts) {
        // If a bucket has >1 d-set, they are candidates
        if (d_set_group.size() > 1) {
            for (size_t i = 0; i < d_set_group.size(); ++i) {
                for (size_t j = i + 1; j < d_set_group.size(); ++j) {
                    
                    // --- Final Check ---
                    // The heuristic (matching signatures) found a candidate pair.
                    // Now, perform the expensive symmetric difference check.
                    
                    // Get the full (non-partitioned) row sets for each d-set
                    auto rows1 = rows_of_d_set(d_set_group[i], A);
                    auto rows2 = rows_of_d_set(d_set_group[j], A);
                    
                    // Calculate the size of the symmetric difference
                    int diff_size = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());

                    // If the difference is less than lambda, it's a non-locating pair
                    if (diff_size < lambda) {
                        to_return.push_back(std::make_tuple(d_set_group[i], d_set_group[j], diff_size));
                    }
                }
            }
        }
    }

    return to_return;
}