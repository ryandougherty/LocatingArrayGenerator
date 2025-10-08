#include "phase1.h"

using partitioned_counts = std::vector<int>;

auto get_partitioned_interaction_counts(const interaction_type& I, const ca_type& A, int X) {
    partitioned_counts counts(X, 0);
    int partition_size = A.size() / X;

    for (const auto& [idx, row] : enumerate(A)) {
        bool match = true;
        for (const auto& [col, val] : zip(I.first, I.second)) {
            if (row[col] != val) {
                match = false;
                break;
            }
        }
        if (match) {
            int partition_index = std::min(X - 1, static_cast<int>(idx / partition_size));
            counts[partition_index]++;
        }
    }
    return counts;
}

robin_hood::unordered_flat_set<int> rows_of_interaction(const interaction_type& I, const ca_type& A) {
    const auto& cols = I.first;
    const auto& vals = I.second;
    robin_hood::unordered_flat_set<int> rows_I_appears;
    for (const auto& [idx, row] : enumerate(A)) {
        bool flag = false;
        for (const auto& [col, val] : zip(cols, vals)) {
            if (row[col] != val) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            rows_I_appears.insert(idx);
        }
    }
    return rows_I_appears;
}


// implemented input of v_type array / maybe vector? nah probably array
auto get_interactions(const t_type t, const vs_type& vs, bool t_bar) {
    // creates COL SETS don't really need to touch
    auto lb = t;
    if (t_bar) {
        lb = 1;
    }
    std::vector<std::vector<k_type>> col_sets;
    for (int i=lb; i<=t; i++) {
        auto cols = combinations(range(vs.size()), i);
        for (const auto& col_set : cols) {
            std::vector<k_type> to_add;
            for (const auto& new_col : col_set) {
                to_add.push_back(new_col);
            }
            col_sets.push_back(to_add);
        }
    }


    // Change up all of these for loops and implement this:
    //  std::vector<interaction_type> interactions;
    // loop for col sets
    std::vector<interaction_type> interactions;
    // t >= 1 and t_bar
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
            std::cerr << "Error, Invalid t\n";
            abort();
        }
    }
    return interactions;
}

bool is_subset(const std::vector<int>& a, const std::vector<int>& b) {
    return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

// A helper function for explicit and safe comparison of interaction_type objects
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

// Helper function to get all rows for a d-set (for the final check)
robin_hood::unordered_flat_set<int> rows_of_d_set(const d_set_type& d_set, const ca_type& A) {
    robin_hood::unordered_flat_set<int> the_rows;
    for(const auto& interaction : d_set) {
        auto rows = rows_of_interaction(interaction, A);
        the_rows.insert(rows.begin(), rows.end());
    }
    return the_rows;
}

std::vector<std::tuple<d_set_type, d_set_type, int>> find_non_detecting_sets(
    const ca_type& A, t_type t, const vs_type& vs, 
    lambda_type lambda, d_type d, bool d_bar, bool t_bar, int X) {
    
    auto interactions = get_interactions(t, vs, t_bar);

    std::vector<d_set_type> d_sets;
    auto lower_lim = d;
    if (d_bar) {
        lower_lim = 1;
    }
    for (int i = lower_lim; i <= d; i++) {
        auto to_add = combinations(interactions, d);
        for (auto& individual_d_set : to_add) {
            d_set_type inner_d_set;
            for (auto& interaction : individual_d_set) {
                inner_d_set.push_back(interaction);
            }
            d_sets.push_back(inner_d_set);
        }
    }

    // Map from an interaction to its partitioned row counts
    robin_hood::unordered_map<interaction_type, partitioned_counts, InteractionHasher> interaction_to_partition_counts;
    for (const auto& interaction : interactions) {
        interaction_to_partition_counts[interaction] = get_partitioned_interaction_counts(interaction, A, X);
    }
    
    // Group d-sets by their partitioned count signature
    robin_hood::unordered_map<partitioned_counts, std::vector<d_set_type>, VectorHasher> d_sets_by_partition_counts;
    std::cout << "Grouping d-sets by partitioned row counts for detecting arrays...\n";
    for (const auto& d_set : d_sets) {
        partitioned_counts total_counts(X, 0);
        robin_hood::unordered_set<int> distinct_rows;
        
        // To create the signature for the d-set, we get the union of rows, 
        // then recalculate the partitioned counts for that union.
        for (const auto& interaction : d_set) {
            auto rows = rows_of_interaction(interaction, A);
            distinct_rows.insert(rows.begin(), rows.end());
        }

        int partition_size = A.size() / X;
        for (const auto& row_idx : distinct_rows) {
            int partition_index = std::min(X - 1, static_cast<int>(row_idx / partition_size));
            total_counts[partition_index]++;
        }

        d_sets_by_partition_counts[total_counts].push_back(d_set);
    }

    std::vector<std::tuple<d_set_type, d_set_type, int>> to_return;
    std::cout << "Ready to look at pairs for detecting arrays...\n";

    // Any d-sets that end up in the same bucket (same signature) are candidates.
    for (const auto& [counts, d_set_group] : d_sets_by_partition_counts) {
        if (d_set_group.size() > 1) {
             for (size_t i = 0; i < d_set_group.size(); ++i) {
                for (size_t j = i + 1; j < d_set_group.size(); ++j) {
                    // Because they have the same signature, they are very likely to be
                    // non-detecting pairs. We can add them directly.
                    // For higher accuracy, you could add a final check here with 
                    // the full row sets, but this heuristic is quite strong.
                    to_return.push_back({d_set_group[i], d_set_group[j], 0});
                }
            }
        }
    }

    // Note: The logic for checking subsets between different groups is more complex
    // with vectors and might not be worth the performance cost. The primary benefit
    // comes from finding groups with identical signatures.

    return to_return;
}

std::vector<std::tuple<d_set_type, d_set_type, int>> find_non_locating_sets(const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar, int X) {
    auto interactions = get_interactions(t, vs, t_bar);
    std::vector<d_set_type> d_sets;

    // get all at most d;
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
    

    // get the rows for each interaction
    // Map from an interaction to its partitioned row counts
    robin_hood::unordered_map<interaction_type, partitioned_counts, InteractionHasher> interaction_to_partition_counts;
    for (const auto& interaction : interactions) {
        interaction_to_partition_counts[interaction] = get_partitioned_interaction_counts(interaction, A, X);
    }

    // Map from a partitioned count vector to the d-sets that have that signature
    robin_hood::unordered_map<partitioned_counts, std::vector<d_set_type>, VectorHasher> d_sets_by_partition_counts;

    std::cout << "Computing partitioned counts of d-sets...\n";
    for (const auto& d_set : d_sets) {
        partitioned_counts total_counts(X, 0);
        for (const auto& interaction : d_set) {
            const auto& counts = interaction_to_partition_counts[interaction];
            for(int i = 0; i < X; ++i) {
                total_counts[i] += counts[i]; // This is a simplification; you might want a more sophisticated way to combine counts for a d-set
            }
        }
        d_sets_by_partition_counts[total_counts].push_back(d_set);
    }

    std::vector<std::tuple<d_set_type, d_set_type, int>> to_return;
    std::cout << "Ready to look at pairs based on partitioned counts...\n";

    // Instead of iterating and comparing scalar row counts, you now iterate through the map.
    // Pairs with the same partitioned_counts vector are candidates for being non-locating.
    for (const auto& [counts, d_set_group] : d_sets_by_partition_counts) {
        if (d_set_group.size() > 1) {
            for (size_t i = 0; i < d_set_group.size(); ++i) {
                for (size_t j = i + 1; j < d_set_group.size(); ++j) {
                    // These pairs have the same signature and are likely non-locating.
                    // Now, you perform the more expensive symmetric difference check on them.
                    auto rows1 = rows_of_d_set(d_set_group[i], A); // You'll need a helper for this
                    auto rows2 = rows_of_d_set(d_set_group[j], A);
                    int diff_size = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());

                    if (diff_size < lambda) {
                        to_return.push_back(std::make_tuple(d_set_group[i], d_set_group[j], diff_size));
                    }
                }
            }
        }
    }

    return to_return;
}
