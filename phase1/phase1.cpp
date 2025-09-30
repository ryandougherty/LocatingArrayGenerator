#include "phase1.h"

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

// template <class InputIterator1, class InputIterator2>
// int size_of_symmetric_difference(InputIterator1 first1, InputIterator1 last1,
//     InputIterator2 first2, InputIterator2 last2)
// {
//     int size = 0;
//     while (true)
//     {
//         if (first1 == last1) return std::distance(first2, last2) + size;
//         if (first2 == last2) return std::distance(first1, last1) + size;

//         if (*first1 < *first2) { ++first1; ++size; }
//         else if (*first2 < *first1) { ++first2; ++size; }
//         else { ++first1; ++first2; }
//     }
// }

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

std::vector<std::tuple<d_set_type, d_set_type, int>> find_non_detecting_sets(const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar) {
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

    robin_hood::unordered_map<interaction_type, robin_hood::unordered_set<int>, InteractionHasher> interaction_to_row_map;
    for (const auto& interaction : interactions) {
        auto rows = rows_of_interaction(interaction, A);
        interaction_to_row_map[interaction] = rows;
    }

    std::map<int, std::vector<std::pair<d_set_type, std::vector<int>>>> d_sets_by_row_count;

    std::cout << "Grouping d-sets by row count...\n";
    for (const auto& d_set : d_sets) {
        robin_hood::unordered_set<int> the_rows;
        for (const auto& interaction : d_set) {
            const auto& rows = interaction_to_row_map[interaction];
            the_rows.insert(rows.begin(), rows.end());
        }
        int n = the_rows.size();
        std::vector<int> vrows(the_rows.begin(), the_rows.end());
        std::sort(vrows.begin(), vrows.end());
        d_sets_by_row_count[n].push_back({d_set, vrows});
    }

    std::vector<std::tuple<d_set_type, d_set_type, int>> to_return;
    std::cout << "Ready to look at pairs for detecting arrays...\n";

    for (auto it1 = d_sets_by_row_count.begin(); it1 != d_sets_by_row_count.end(); ++it1) {
        for (auto it2 = it1; it2 != d_sets_by_row_count.end(); ++it2) {
            const auto& group1 = it1->second;
            const auto& group2 = it2->second;

            if (it1 == it2) { // Same size group
                for (size_t i = 0; i < group1.size(); ++i) {
                    for (size_t j = i + 1; j < group1.size(); ++j) {
                        if (group1[i].second == group1[j].second) { // Identical row sets
                            auto d_set1 = group1[i].first;
                            auto d_set2 = group1[j].first;

                            // Sort using the new, explicit comparison function
                            std::sort(d_set1.begin(), d_set1.end(), compare_interactions);
                            std::sort(d_set2.begin(), d_set2.end(), compare_interactions);

                            // Only add the pair if they are not permutations of each other
                            if (d_set1 != d_set2) {
                                to_return.push_back({group1[i].first, group1[j].first, 0});
                            }
                        }
                    }
                }
            }
            else { // Different size groups
                for (const auto& pair1 : group1) {
                    for (const auto& pair2 : group2) {
                        // Check if the smaller is a subset of the larger
                        if (is_subset(pair1.second, pair2.second)) {
                            to_return.push_back({pair1.first, pair2.first, 0});
                        }
                    }
                }
            }
        }
    }
    return to_return;
}

std::vector<std::tuple<d_set_type, d_set_type, int>> find_non_locating_sets(const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar) {
    auto interactions = get_interactions(t, vs, t_bar);

    // get all at most d;
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
    

    // get the rows for each interaction
    robin_hood::unordered_map<interaction_type, robin_hood::unordered_set<int>, InteractionHasher> interaction_to_row_map;
    for (const auto& interaction : interactions) {
        auto rows = rows_of_interaction(interaction, A);
        interaction_to_row_map[interaction] = rows;
    }


    
    // get the (last) lambda rows for each d_set. Will check below all pairs from sets that have less than lambda in common (if not, then they must have lambda or more in symmetric difference).

    std::map<int, std::vector<std::pair<d_set_type, std::vector<int>>>> initial_rows_map;

    std::cout << "Computing rows of d-sets...\n";
    
    std::vector<int> dest(lambda);
    for (const auto& d_set : d_sets) {
        d_set_type copied_d_set;
        for (const auto& interaction : d_set) {
            copied_d_set.push_back(interaction);
        }
        robin_hood::unordered_set<int> the_rows;
        for (const auto& interaction : copied_d_set) {
            const auto& rows = interaction_to_row_map[interaction];

            the_rows.insert(rows.begin(), rows.end());
        }
        int n = the_rows.size();
        std::vector<int> vrows(the_rows.begin(), the_rows.end());
        std::sort(vrows.begin(), vrows.end());

        if (initial_rows_map.count(n)) {
            initial_rows_map[n].push_back(std::make_pair(copied_d_set, vrows));
        }
        else {
            std::vector<std::pair<d_set_type, std::vector<int>>> the_inner_vector{std::make_pair(copied_d_set, vrows)};
            auto s = std::make_pair(n, the_inner_vector);
            initial_rows_map.insert(s);
        }
    }
    std::vector<std::pair<int, std::vector<std::pair<d_set_type, std::vector<int>>>>> largest_rows_num_map(initial_rows_map.begin(), initial_rows_map.end());
    std::sort(largest_rows_num_map.begin(), largest_rows_num_map.end());
    for (auto& [key, vec_of_inner_pairs] : largest_rows_num_map) {
        std::sort(vec_of_inner_pairs.begin(), vec_of_inner_pairs.end(), 
            [](auto& pair1, auto& pair2) { return pair1.second < pair2.second;  }
        );
    }
  
    // SORTED vector for the first AND second values
    // std::vector<std::pair<std::vector<int>, std::vector<std::pair<d_set_type, std::vector<int>>>>> largest_rows_num_map;
    // 
    // using interaction_type = std::pair<std::vector<int>, std::vector<int>>;
    // using d_set_type = std::vector<interaction_type>;
    // 
    // 
    // iterate through all pairs of row_nums, and only consider those that have less than lambda symm diff

    // the inner tuple is the d-set pair, and how many rows THEY HAVE ALREADY BEEN SEPARATED
    std::vector<std::tuple<d_set_type, d_set_type, int>> to_return;
    std::cout << "Ready to look at pairs...\n";
    std::cout << "largest_rows_num_map size=" << largest_rows_num_map.size() << "\n";

    for (const auto& pair : largest_rows_num_map) {
        std::cout << "(" << pair.first << ", " << pair.second.size() << ") ";
    }
    std::cout << "\n";

    // largest_rows_num_map is SORTED
    //      so if the nums are sufficiently far apart, all pairs of them must be locating
    //      i.e., pair2.num_rows - pair1.num_rows >= lambda
    for (auto pair1 = largest_rows_num_map.begin(); pair1 != largest_rows_num_map.end(); pair1++) {
        const auto& num_rows1 = (*pair1).first;
        const auto& all_dset1 = (*pair1).second;

        for (auto pair2 = pair1; pair2 != largest_rows_num_map.end(); pair2++) {
            const auto& num_rows2 = (*pair2).first;
            const auto& all_dset2 = (*pair2).second;

            // int symm_size = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());
            if (num_rows2 - num_rows1 >= lambda) {
                break;
            }

            if (num_rows1 == num_rows2) {
                // go over all pairs without repeats
                for (auto inner_pair1 = all_dset1.begin(); inner_pair1 != all_dset1.end(); inner_pair1++) {
                    const auto& [d_set1, rows_1] = *inner_pair1;

                    for (auto inner_pair2 = std::next(inner_pair1); inner_pair2 != all_dset1.end(); inner_pair2++) {
                        const auto& [d_set2, rows_2] = *inner_pair2;

                        if (rows_1.back() < rows_2.front()) {
                            break;
                        }

                        if (d_set1 == d_set2) {
                            continue;
                        }

                        int diff_size = size_of_symmetric_difference(rows_1.begin(), rows_1.end(), rows_2.begin(), rows_2.end());

                        if (diff_size < lambda) {
                            to_return.push_back(std::make_tuple(d_set1, d_set2, diff_size));
                        }
                    }
                }
            } else {
                for (const auto& [d_set1, rows_1] : all_dset1) {
                    for (const auto& [d_set2, rows_2] : all_dset2) {
                        if (rows_1.back() < rows_2.front()) {
                            break;
                        }
                        int diff_size = size_of_symmetric_difference(rows_1.begin(), rows_1.end(), rows_2.begin(), rows_2.end());
                        if (diff_size < lambda) {
                            to_return.push_back(std::make_tuple(d_set1, d_set2, diff_size));
                        }
                    }
                }
            }

        }
    }

    return to_return;
}