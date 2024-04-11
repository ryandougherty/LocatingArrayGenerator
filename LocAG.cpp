#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <map>
#include <random>
#include <ranges>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>


#include "combinations.hpp"
#include "enumerate.hpp"
//#include "flat_hash_map.hpp"
#include "product.hpp"
#include "range.hpp"
#include "robin_hood.h"
#include "zip.hpp"
#include "LocAG.h"

using N_type = int;
using d_type = uint8_t;
using k_type = uint8_t;
using v_type = uint8_t;
using vs_type = std::vector<v_type>;
using t_type = uint8_t;
using lambda_type = uint8_t;

using ca_type = std::vector<std::vector<v_type>>;
using interaction_type = std::pair<std::vector<k_type>, std::vector<v_type>>;
using d_set_type = std::vector<interaction_type>;

using ga_individual_type = std::vector<double>;
using ga_fitness_type = std::tuple<ga_individual_type, int, long long>;

using namespace iter;
using namespace std::chrono;

std::mt19937_64 rng(0); // set seed
std::uniform_real_distribution<double> unif(0, 1);
std::uniform_int_distribution<int> ind_size(10, 30);
std::uniform_int_distribution<int> any_int;

auto interaction_to_str(const interaction_type& I) {
    std::string result = "(";
    for (auto& col : I.first) {
        result += std::to_string(col) + ',';
    }
    result += "),(";
    for (auto& val : I.second) {
        result += std::to_string(val) + ',';
    }
    return result + ")";
}

void print_interaction(const interaction_type& I) {
    std::cout << interaction_to_str(I);
}

void print_d_set(const d_set_type& D) {
    for (const auto& I : D) {
        print_interaction(I);
    }
}

auto d_set_to_str(const d_set_type& D) {
    std::string result;
    for (const auto& I : D) {
        result += interaction_to_str(I) + " ";
    }
    return result;
}

void print_array(const ca_type& A) {
    for (const auto& row : A) {
        for (const auto& value : row) {
            std::cout << std::to_string(value) << ' ';
        }
        std::cout << '\n';
    }
}

template <typename T>
void print_vec(std::vector<T> I) {
    for (const auto& elem : I) {
        std::cout << std::to_string(elem) << ",";
    }
}


long long comb(unsigned n, unsigned k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    long long result = n;
    for (int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

double calc_p(const int t, const int v) {
    return 1 / pow(v, t);
}

ca_type random_array(const N_type N, const k_type k, const vs_type& vs) {
    ca_type to_return(N, std::vector<v_type>(k, 0));
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < k; col++) {
            std::uniform_int_distribution<int> distribution(0, vs[col] - 1);
            to_return[row][col] = distribution(rng);
        }
    }
    return to_return;
}

// Source: https://stackoverflow.com/questions/10405030/c-unordered-map-fail-when-used-with-a-vector-as-key
struct VectorHasher {
    int operator()(const std::vector<v_type>& V) const {
        int hash = V.size();
        for (auto& i : V) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};
struct InteractionHasher {
    int operator()(const interaction_type& V) const {
        int hash = V.first.size();
        for (auto& i : V.first) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        for (auto& i : V.second) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};
struct DSetHasher {
    int operator()(const d_set_type& V) const {
        int hash = V.size();
        for (auto& I : V) {
            for (auto& i : I.first) {
                hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            for (auto& i : I.second) {
                hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
        }
        
        return hash;
    }
};

auto first_uncovered_cols(ca_type A, const t_type t, const k_type k, const vs_type& vs, const lambda_type lambda) {
    std::vector<v_type> row_in_A(t, 0);
    std::vector<k_type> cols_to_return;
    for (const auto& cols : combinations(range(k), t)) {
        robin_hood::unordered_flat_map<std::vector<v_type>, int, VectorHasher> c;
        for (const auto& row : A) {
            for (int i = 0; i < cols.size(); i++) {
                row_in_A[i] = row[cols[i]];
            }
            if (c.find(row_in_A) != c.end()) {
                c[row_in_A] += 1;
            }
            else {
                c[row_in_A] = 1;
            }
        }
        if (c.size() != pow(vs[cols.size()], t)) {
            for (const auto& col : cols) {
                cols_to_return.push_back(col);
            }
            return cols_to_return;
        }
        else {
            const auto& it = *std::min_element(std::begin(c), std::end(c),
                [](const auto& l, const auto& r) { return l.second < r.second; });
            if (it.second < lambda) {
                for (const auto& col : cols) {
                    cols_to_return.push_back(col);
                }
                return cols_to_return;
            }
        }
    }
    return std::vector<k_type>();
}

auto rows_of_interaction(const interaction_type& I, const ca_type& A) {
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

template <class InputIterator1, class InputIterator2>
int size_of_symmetric_difference(InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2)
{
    int size = 0;
    while (true)
    {
        if (first1 == last1) return std::distance(first2, last2) + size;
        if (first2 == last2) return std::distance(first1, last1) + size;

        if (*first1 < *first2) { ++first1; ++size; }
        else if (*first2 < *first1) { ++first2; ++size; }
        else { ++first1; ++first2; }
    }
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

auto find_non_locating_sets(const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar) {
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

auto read_ca_from_cagen(const std::string& filename, const vs_type& vs) {
    std::ifstream file(filename);
    std::string line;
    // skip first line as that does not have CA rows
    std::getline(file, line);
    ca_type result;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<std::string> things_in_line;
        while (ss.good()) {
            std::string substr;
            std::getline(ss, substr, ',');
            things_in_line.push_back(substr);
        }
        std::vector<v_type> ca_line;
        for (const auto& [col_idx, elem] : enumerate(things_in_line)) {
            if (elem == "*") {
                ca_line.push_back(any_int(rng) % vs[col_idx]);
            } else {
                ca_line.push_back(std::stoi(elem));
            }
        }

        result.push_back(ca_line);
    }
    return result;
}

auto lookup_or_assign_interaction_map(std::map<interaction_type, robin_hood::unordered_flat_set<int>>& rows_map, const interaction_type& interaction, const ca_type& ind) {
    if (rows_map.find(interaction) != rows_map.end()) {
        const auto& rows = rows_map[interaction];
        return rows;
    } else {
        const auto& rows = rows_of_interaction(interaction, ind);
        rows_map[interaction] = rows;
        return rows;
    } 
}

int fitness(const ca_type& ind, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, const int& threshold) {
    int score = 0;
    // value is sorted
    std::unordered_map<d_set_type, std::vector<N_type>, DSetHasher> rows_map;
    auto rows_of_dset = [=,&rows_map](const d_set_type& d_set) {
        if (rows_map.find(d_set) != rows_map.end()) {
            return rows_map[d_set];
        } else {
            robin_hood::unordered_set<N_type> the_rows;
            for (const auto& interaction : d_set) {
                const auto& rows = rows_of_interaction(interaction,ind);
                the_rows.insert(rows.begin(), rows.end());
            }
            std::vector<int> vrows(the_rows.begin(), the_rows.end());
            std::sort(vrows.begin(), vrows.end());
            rows_map[d_set] = vrows;
            return vrows;
        }
    };
    for (const auto& [dset_1, dset_2, num_times_sep_already] : non_locating_pairs) {
        auto requirement = l - num_times_sep_already;
        auto rows1 = rows_of_dset(dset_1);
        auto rows2 = rows_of_dset(dset_2);
        int n = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());
      
        if (n >= requirement) {
            score += 1;
        }

        if (score >= threshold) {
            return threshold+1;
        }
    }
    return score;
}

ca_type cross(const ca_type& p1, const ca_type& p2, d_type d, t_type t, const vs_type& vs, lambda_type l) {
    int val = any_int(rng) % 2;
    int n = p1.size();
    ca_type child;
    if (val == 0) {
        auto rand_idx = any_int(rng) % p1.size();
        for (int i=0; i<rand_idx; i++) {
            child.push_back(p1[i]);
        }
        for (int i=rand_idx; i<n; i++) {
            child.push_back(p2[i]);
        }
    } else if (p1.size() != 1) {
        auto rand_idx1 = any_int(rng) % (p1.size());
        auto rand_idx2 = any_int(rng) % (p1.size());
        while (rand_idx1 == rand_idx2) {
            rand_idx1 = any_int(rng) % (p1.size());
            rand_idx2 = any_int(rng) % (p1.size());
        }
        auto lower = std::min(rand_idx1,rand_idx2);
        auto higher = std::max(rand_idx1,rand_idx2);
        for (int i=0; i<lower; i++) {
            child.push_back(p1[i]);
        }
        for (int i=lower; i<higher; i++) {
            child.push_back(p2[i]);
        }
        for (int i=higher; i<p1.size(); i++) {
            child.push_back(p1[i]);
        }
    } else {
        val = any_int(rng) % 2;
        if (val == 0) {
            child = p1;
        } else {
            child = p2;
        }
    }
    return child;
}

ca_type mutate(const ca_type& p1, d_type d, t_type t, const vs_type& vs, lambda_type l) {
    int val = any_int(rng) % 3;
    int n = p1.size();
    ca_type child = p1;
    if (val == 0) {
        auto rand_row = any_int(rng) % n;
        for (int col=0; col<p1[0].size(); col++) {
            child[rand_row][col] = any_int(rng) % vs[col];
        }
    } else if (val == 1) {
        auto rand_col = any_int(rng) % vs.size();
        for (int row=0; row<n; row++) {
            child[row][rand_col] = any_int(rng) % vs[rand_col];
        }
    } else {
        auto rand_row = any_int(rng) % p1.size();
        auto rand_col = any_int(rng) % p1[0].size();
        auto rand_val = any_int(rng) % vs[rand_col];
        child[rand_row][rand_col] = rand_val;
    }
    return child;
}

struct Ind_NonRecompute_Fitness {
    ca_type A;
    int fitness;
};

struct PercentGAFitnessInd {
    std::vector<double> percents;
    int N = -1;
    long long time = -1;

    bool operator==(PercentGAFitnessInd const&) const = default;
};

ca_type try_N(N_type N, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, double percent) {

    ca_type s;
    int pop_size = 100;
    int num_gens = 50;

    std::vector<Ind_NonRecompute_Fitness> pop(pop_size);
     
    for (auto& elem : pop) {
        elem.A = random_array(N, vs.size(), vs);
        elem.fitness = -1;
    }

    auto best_fitness = std::numeric_limits<int>::min();
    const int max_possible_fitness = non_locating_pairs.size() * percent;

    for (int gen=0; gen<num_gens; gen++) {
        std::vector<std::pair<int, Ind_NonRecompute_Fitness>> fitnesses;
        for (auto& I : pop) {
            int f = I.fitness;
            if (f == -1) {
                f = fitness(I.A, d, t, vs, l, non_locating_pairs, max_possible_fitness);
                I.fitness = f; 
            }
            if (f >= max_possible_fitness) {
                return I.A;
            }
            fitnesses.push_back(std::make_pair(f,I));
        }


        std::sort(fitnesses.begin(), fitnesses.end(), [](const auto& first, const auto& second) {
            return first.first < second.first;
        });

        std::vector<Ind_NonRecompute_Fitness> new_vec;
        for (int i=pop_size/2; i < pop_size; i++) {
            new_vec.push_back(fitnesses[i].second);
        }
        pop = new_vec;
        new_vec.clear();

        while (new_vec.size() < pop_size / 2) {
            auto idx1 = any_int(rng) % pop.size();
            auto idx2 = any_int(rng) % pop.size();

            const auto& p1 = pop[idx1];
            const auto& p2 = pop[idx2];

            auto cross_percent = any_int(rng) % 10;
            auto mut_percent = any_int(rng) % 10;
            if (cross_percent == 0 && mut_percent < 3) {
                auto new_ind = cross(p1.A,p2.A,d,t,vs,l);
                new_ind = mutate(new_ind,d,t,vs,l);
                Ind_NonRecompute_Fitness true_new_ind;
                true_new_ind.A = new_ind;
                true_new_ind.fitness = -1;
                new_vec.push_back(true_new_ind);
            } else if (cross_percent == 0) {
                auto new_ind = cross(p1.A,p2.A,d,t,vs,l);
                Ind_NonRecompute_Fitness true_new_ind;
                true_new_ind.A = new_ind;
                true_new_ind.fitness = -1;
                new_vec.push_back(true_new_ind);
            }
        }
        pop.insert(pop.end(), new_vec.begin(), new_vec.end());
    }

    return s;
}




ca_type try_N_SA(N_type N, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& only_these_pairs) {



    ca_type empty;
    ca_type A = random_array(N, vs.size(), vs);
    auto temp = 1.0;
    auto rate = 0.99;
    auto num_iter = 1000;

    auto required_fitness = only_these_pairs.size();

    auto f = fitness(A, d, t, vs, l, only_these_pairs, required_fitness);
    for (int it=0; it<num_iter; it++) {
        
        if (f >= required_fitness) {
            return A;
        }
        auto A_prime = mutate(A, d, t, vs, l);
        auto f_prime = fitness(A_prime, d, t, vs, l, only_these_pairs, required_fitness);
        // std::cout << "req=" << required_fitness << ", got=" << f_prime << "\n";
        auto diff = f_prime - f;
        if (f_prime >= f) {
            A = A_prime;
            f = f_prime;
        } else {
            auto prob = std::exp(-diff / temp);
            if (prob < unif(rng)) {
                A = A_prime;
                f = f_prime;
            }
        }
        temp = rate * temp;
    }

    if (f >= required_fitness) {
        return A;
    }

    // std::cout << "Needed " << required_fitness << ", got " << f << "\n";
    // print_array(A);
    return empty;
}





// insert new parameter
// parameter for the percentage of completion of locating rows

ca_type go(const d_type& d, const t_type& t, const vs_type& vs, const lambda_type& l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, const double& percent) {
    
    
    bool succ_first = true;
    ca_type result;
    std::size_t howManyPairs;

    const std::vector<std::tuple<d_set_type, d_set_type, int>> only_these_pairs(non_locating_pairs.begin(), non_locating_pairs.begin() + non_locating_pairs.size() * percent);
    
    int N = 1;
    for (const auto& [d_set1, d_set2, num] : only_these_pairs) {
        N = std::max(N, l-num);
    }

    while (true) {

        result = try_N_SA(N, d, t, vs, l, only_these_pairs);
        if (succ_first &&  result.size() > 0) {
            return result;
        }
        if (result.size() > 0) {
            break;
        }
        N *= 2;
        // std::cout << "Trying " << N << "\n";
        succ_first = false;
    }

    int N_hi = N;
    int N_lo = N / 2;
    while (N_lo < N_hi) {
        int N_mid = (N_lo + N_hi) / 2;
        auto result2 = try_N_SA(N_mid, d, t, vs, l, only_these_pairs);
        if (result2.size() > 0) {
            N_hi = N_mid;
            result = result2;
        } else {
            N_lo = N_mid + 1;
        }
    }
    return result;
}

auto parse_vs(const std::vector<std::string>& exp_params) {
    vs_type new_params;
    for (const auto& str : exp_params) {
        auto idx = str.find("^");
        v_type v = std::atoi(std::string(str.begin(), str.begin() + idx).c_str());
        auto k = std::atoi(std::string(str.begin() + idx+1, str.end()).c_str());
        for (int i=0; i<k; i++) {
            new_params.push_back(v);
        }
    }
    return new_params;
}




bool dominates(const PercentGAFitnessInd& ind, const PercentGAFitnessInd& other) {
    return (ind.N <= other.N && 
        ind.time <= other.time);
}

auto pareto_and_rest(std::vector<PercentGAFitnessInd> points) {
    int candidate_ind_number = 0;
    std::vector<PercentGAFitnessInd> dominated_pts;
    std::vector<PercentGAFitnessInd> pareto;
    while (true) {
        auto candidate_ind = points[candidate_ind_number];
        points.erase(points.begin() + candidate_ind_number);
        bool non_dominated = true; 
        int ind_number = 0;
        while (points.size() != 0 && ind_number < points.size()) {
            auto ind = points[ind_number];
            if (dominates(candidate_ind, ind)) {
                points.erase(points.begin() + ind_number);
                dominated_pts.push_back(ind);
            } else if (dominates(ind, candidate_ind)) {
                non_dominated = false;
                dominated_pts.push_back(candidate_ind);
                ind_number++;
            } else {
                ind_number++;
            }
        }
        if (non_dominated) {
            pareto.push_back(candidate_ind);
        }
        if (points.size() == 0) {
            break;
        }
    }
    return std::make_pair(pareto, dominated_pts);
}





auto generate_rand_percent_individual() {
    std::vector<double> percents;
    int rand_length = ind_size(rng); // between 10 and 30 stages
    for (int i=0; i<rand_length; i++) {
        percents.push_back(unif(rng));
    }
    std::sort(percents.begin(), percents.end());
    percents[0] = 0.001;
    percents.push_back(1.0);
    PercentGAFitnessInd result;
    result.percents = percents;
    return result;
}



auto percent_GA(d_type d, t_type t, const vs_type& vs, const lambda_type& l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, bool use_default_percents) {

    if (use_default_percents) {
        const std::vector<double> percents = {0.021576,0.021576,0.022644,0.030792,0.090424,0.071014,0.083679,0.172455,0.220123,0.415283,1.000000};
//{0.001000,0.002625,0.010040,0.017456,0.631592,0.094666,0.151855,0.167999,0.172241,1.000000};
        auto non_locating_pairs_copy = non_locating_pairs;
        int num_rows = 0;
        auto start = high_resolution_clock::now();
        for (const auto& percent : percents) { 
            std::vector<std::tuple<d_set_type, d_set_type, int>> new_non_locating_pairs;

            auto ga_rows = go(d,t,vs,l,non_locating_pairs_copy,percent);
            num_rows += ga_rows.size();

            for (const auto& [dset_1, dset_2, num_times_sep_already] : non_locating_pairs_copy) {
                std::unordered_map<d_set_type, std::vector<N_type>, DSetHasher> rows_map;
                auto rows_of_dset = [=,&rows_map](const d_set_type& d_set) {
                    if (rows_map.find(d_set) != rows_map.end()) {
                        return rows_map[d_set];
                    } else {
                        robin_hood::unordered_set<N_type> the_rows;
                        for (const auto& interaction : d_set) {
                            const auto& rows = rows_of_interaction(interaction,ga_rows);
                            the_rows.insert(rows.begin(), rows.end());
                        }
                        std::vector<int> vrows(the_rows.begin(), the_rows.end());
                        std::sort(vrows.begin(), vrows.end());
                        rows_map[d_set] = vrows;
                        return vrows;
                    }
                }; 

                auto rows1 = rows_of_dset(dset_1);
                auto rows2 = rows_of_dset(dset_2);
                int n = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());
                if (num_times_sep_already + n < l) {
                    new_non_locating_pairs.push_back(std::make_tuple(dset_1, dset_2, num_times_sep_already + n));
                }
            }

            non_locating_pairs_copy = new_non_locating_pairs;

            if (non_locating_pairs_copy.size() == 0) {
                break;
            }

            std::cout << "Added " << num_rows << "rows, there are " << non_locating_pairs_copy.size() << " remaining pairs\n";
            new_non_locating_pairs.clear();
        }
        auto stop = high_resolution_clock::now();
        auto total_time = duration_cast<milliseconds>(stop-start).count();
        std::vector<PercentGAFitnessInd> result;
        PercentGAFitnessInd ind;
        ind.N = num_rows;
        ind.percents = percents;
        ind.time = total_time;
        result.push_back(ind);
        return result;
    }

    int pop_size = 100;
    int num_gens = 50;

    std::vector<PercentGAFitnessInd> result;

    std::vector<PercentGAFitnessInd> pop;
    for (int i=0; i<pop_size; i++) {
        pop.push_back(generate_rand_percent_individual());
    }

    for (int gen=0; gen<num_gens; gen++) {

        std::cout << "Generation #" << gen << "\n";

        std::vector<PercentGAFitnessInd> fitnesses;
        int individual = 0;
        for (auto& I : pop) {
            // std::cout << "Fitness for individual #" << individual << ": ";
            

            int num_rows = 0;
            long long ind_total_time = 0;
            if (I.N != -1 && I.time != -1) {
                num_rows = I.N;
                ind_total_time = I.time;
            } else {
                auto percents = I.percents;
                auto non_locating_pairs_copy = non_locating_pairs;

                auto start = high_resolution_clock::now();
                for (const auto& percent : percents) { 

                    std::vector<std::tuple<d_set_type, d_set_type, int>> new_non_locating_pairs;

                    auto ga_rows = go(d,t,vs,l,non_locating_pairs_copy,percent);
                    num_rows += ga_rows.size();

                    for (const auto& [dset_1, dset_2, num_times_sep_already] : non_locating_pairs_copy) {
                        std::unordered_map<d_set_type, std::vector<N_type>, DSetHasher> rows_map;
                        auto rows_of_dset = [=,&rows_map](const d_set_type& d_set) {
                            if (rows_map.find(d_set) != rows_map.end()) {
                                return rows_map[d_set];
                            } else {
                                robin_hood::unordered_set<N_type> the_rows;
                                for (const auto& interaction : d_set) {
                                    const auto& rows = rows_of_interaction(interaction,ga_rows);
                                    the_rows.insert(rows.begin(), rows.end());
                                }
                                std::vector<int> vrows(the_rows.begin(), the_rows.end());
                                std::sort(vrows.begin(), vrows.end());
                                rows_map[d_set] = vrows;
                                return vrows;
                            }
                        }; 

                        auto rows1 = rows_of_dset(dset_1);
                        auto rows2 = rows_of_dset(dset_2);
                        int n = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());
                        if (num_times_sep_already + n < l) {
                            new_non_locating_pairs.push_back(std::make_tuple(dset_1, dset_2, num_times_sep_already + n));
                        }
                    }

                    non_locating_pairs_copy = new_non_locating_pairs;
                    new_non_locating_pairs.clear();
                }
                auto stop = high_resolution_clock::now();
                ind_total_time = duration_cast<milliseconds>(stop-start).count();
                I.N = num_rows;
                I.time = ind_total_time;
            }
            // std::cout << individual << "," << num_rows << "," << ind_total_time << ",I=";
            individual++;
            // print_vec(I.percents);
            // std::cout << "\n";
            fitnesses.push_back(I);
        }
        auto target_size = pop_size/2;
        auto [pareto, rest] = pareto_and_rest(fitnesses);
        for (auto& [percents, N, time] : pareto) {
            std::cout << N << "," << time << ",";
            print_vec(percents);
            std::cout << "\n";
        }
        std::vector<PercentGAFitnessInd> new_pop(pareto.begin(), pareto.end());

        result = pareto; // save results each time

        while (new_pop.size() < target_size) {
            for (auto& elem : pareto) {
                fitnesses.erase(std::remove(fitnesses.begin(), fitnesses.end(), elem), fitnesses.end());
            }
            auto [pareto2, rest2] = pareto_and_rest(fitnesses); 
            for (auto& elem : pareto2) {
                if (new_pop.size() < target_size) {
                    new_pop.push_back(elem);
                } else {
                    break;
                }
            }
            pareto = pareto2;
        }
        pop.clear();
        for (auto& elem : new_pop) {
            pop.push_back(elem); // the individual
        }
        new_pop.clear();

        std::vector<PercentGAFitnessInd> children;

        // crossover/mutate
        while (children.size() < pop_size/2) {
            // crossover 
            auto rand_p1 = pop[any_int(rng) % pop.size()];
            auto rand_p2 = pop[any_int(rng) % pop.size()];
            auto rand_idx = any_int(rng) % std::min(rand_p1.percents.size(), rand_p2.percents.size());
            while (rand_idx == 0) {
                rand_idx = any_int(rng) % std::min(rand_p1.percents.size(), rand_p2.percents.size());
            }

            PercentGAFitnessInd child;
            for (int i=0; i<rand_idx; i++) {
                child.percents.push_back(rand_p1.percents[i]);
            }
            for (int i=rand_idx; i<rand_p2.percents.size(); i++) {
                child.percents.push_back(rand_p2.percents[i]);
            }
            children.push_back(child);
        }

        // mutate
        for (auto& child : children) {
            auto r = any_int(rng) % 10;
            if (r == 0) {
                // split
                auto rand_idx = any_int(rng) % (child.percents.size()-1); // ensure not the last index
                auto rand_val = unif(rng);
                child.percents.insert(child.percents.begin() + rand_idx, rand_val);
            } else if (r == 1) {
                // join
                auto rand_idx = any_int(rng) % (child.percents.size()-1); // ensure not the last index 
                child.percents.erase(child.percents.begin() + rand_idx);
            } else if (r == 2) {
                // mutate one entry randomly
                auto rand_idx = any_int(rng) % (child.percents.size()-1); // ensure not the last index
                // again terrible
                child.percents[rand_idx] = unif(rng);
            }

            pop.push_back(child);
        }
    }

    return result;
}


// the reason these are not all ascending order of # levels is that we did this
//          before critical thinking was invented.
const std::unordered_map<std::string, std::vector<std::string>> configs {
    {"Apache", {"2^158", "3^8", "4^4", "5^1", "6^1"}},
    {"Bugzilla", {"2^49", "3^1", "4^2"}},
    {"Flex", {"5^2", "3^4", "2^23"}}, // done
    {"GCC", {"2^189", "3^10"}},
    {"Make", {"6^1", "5^1", "4^2", "3^4", "2^14"}},
    {"Mobile", {"10^8", "9^1", "8^4", "7^5", "6^10", "5^4", "4^6", "3^9", "2^28"}},
    {"SPINS", {"2^13", "4^5"}}, // done
    {"SPINV", {"2^42", "3^2", "4^11"}},
    {"TCAS", {"2^7", "3^2", "4^1", "10^4"}},
    {"Wireless", {"5^9", "4^5", "3^7", "2^3"}}
};

auto lookup_config_and_params(const std::string& config_name, const t_type t, const lambda_type lambda) {

    std::string prefix = "./evaluation/";
    vs_type the_vals;
    if (configs.find(config_name) != configs.end()) {
        the_vals = parse_vs(configs.at(config_name));
        auto filename = prefix + config_name + "/" + config_name + "_" + std::to_string(t) + "_" + std::to_string(lambda) + ".csv";
        return std::make_pair(the_vals, filename);
    } else {
        std::vector<std::string> the_new_config{config_name};
        the_vals = parse_vs(the_new_config);
        auto filename = prefix + config_name + "-t" + std::to_string(t) + "_l" + std::to_string(lambda) + ".csv";
        return std::make_pair(the_vals, filename);
    }
    
}

int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Usage: ./LocAG <name of config>\n";
        return -1;
    }

    

    for (d_type d = 1; d <= 1; d++) {
        for (t_type t = 2; t <= 2; t++) {
            for (lambda_type lambda = 1; lambda <= 4; lambda++) {
                const std::string config_name = argv[1];
                const auto& [vs, filename] = lookup_config_and_params(config_name, t, lambda);

                const bool d_bar = true;
                const bool t_bar = true;

                assert(d < *std::min_element(vs.begin(), vs.end()));

                std::cout << "------------d=" << std::to_string(d) << ", t=" << std::to_string(t) << ", lambda=" << std::to_string(lambda) << ", filename=" << filename << "------------\n";

                auto start = high_resolution_clock::now();
                ca_type A = read_ca_from_cagen(filename, vs);
                // std::cout << "Read file with " << A.size() << " rows.\n";

                // Finds initial non_locating_pairs
                auto non_locating_pairs = find_non_locating_sets(A, t, vs, lambda, d, d_bar, t_bar);
                auto stop = high_resolution_clock::now();
                auto first_stage_N = A.size();
                auto first_stage_time = duration_cast<milliseconds>(stop-start).count();

                std::cout << "First Stage N=" << first_stage_N << ", Time=" << first_stage_time << "\n"; 

                std::cout << "There are " << non_locating_pairs.size() << " remaining non-locating pairs\n";

                /* ------------------------------------- Stage 2: GA ------------------------------ */

                if (non_locating_pairs.size() == 0) {
                    break;
                }

                auto pareto = percent_GA(d,t,vs,lambda,non_locating_pairs,true);
                for (const auto& [percents, num_rows, time] : pareto) {
                    std::cout << "N total=" << first_stage_N + num_rows << ", Time total=" << first_stage_time + time << ", percents=";
                    print_vec(percents);
                    std::cout << "\n";
                }
                
            }
        }
    }
    

}
