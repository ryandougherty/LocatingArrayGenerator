// LocatingDetectingGeneticAlg.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <vector>


#include "combinations.hpp"
#include "enumerate.hpp"
//#include "flat_hash_map.hpp"
#include "product.hpp"
#include "range.hpp"
#include "robin_hood.h"
#include "zip.hpp"
#include "LocatingDetectingGeneticAlg.h"

using N_type = int;
using d_type = uint8_t;
using k_type = uint8_t;
using v_type = uint8_t;
using t_type = uint8_t;
using lambda_type = uint8_t;

using ca_type = std::vector<std::vector<v_type>>;
using interaction_type = std::pair<std::vector<k_type>, std::vector<v_type>>;
using d_set_type = std::vector<interaction_type>;

using namespace iter;
using namespace std::chrono;



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
    //for (int i = 0; i < I.first.size(); i++) {
    //    auto col = I.first[i];
    //    auto val = I.second[i];
    //    std::cout << '(' << col << ', ' << val << '),';
    //}
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

long long LLL_CA_sum(const int t, const int k, const int v, const int lambda) {
    const double p = calc_p(t, v);
    long long N = 1;
    auto calc_quantity = [t, k, v, lambda, p](long long m) {
        const auto e = std::exp(1.0);
        auto sum = 0.0;
        for (int i = 0; i < lambda; i++) {
            sum += comb(m, i) * pow(p, i) * pow(1 - p, m - i);
        }
        return e * (comb(k, t) - comb(k - t, t)) * pow(v, t) * sum;
    };
    while (true) {
        auto q = calc_quantity(N);
        if (q > 1) {
            N *= 2;
        }
        else {
            break;
        }
    }
    auto lb = N / 2;
    auto ub = N;
    while (lb < ub) {
        auto m = std::midpoint(lb, ub);
        if (calc_quantity(m) > 1) {
            lb = m + 1;
        }
        else if (calc_quantity(m + 1) > 1) {
            return m;
        }
        else {
            ub = m;
        }
    }
    return lb;
}

ca_type random_array(const N_type N, const k_type k, const v_type v) {
    std::random_device generator;
    std::uniform_int_distribution<int> distribution(0, v - 1);
    ca_type to_return(N, std::vector<v_type>(k, 0));
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < k; col++) {
            to_return[row][col] = distribution(generator);
        }
    }
    return to_return;
}

// https://stackoverflow.com/questions/10405030/c-unordered-map-fail-when-used-with-a-vector-as-key
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

auto first_uncovered_cols(ca_type A, const t_type t, const k_type k, const v_type v, const lambda_type lambda) {
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
        if (c.size() != pow(v, t)) {
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

ca_type LLL_gen(t_type t, k_type k, v_type v, lambda_type lambda) {
    std::random_device generator;
    std::uniform_int_distribution<int> distribution(0, v - 1);
    auto N = LLL_CA_sum(t, k, v, lambda);
    auto A = random_array(N, k, v);
    int i = 0;
    while (true) {
        i += 1;
        auto cols = first_uncovered_cols(A, t, k, v, lambda);
        if (!cols.empty()) {
            for (auto& row : A) {
                for (auto& col : cols) {
                    row[col] = distribution(generator);
                }
            }
        }
        else {
            return A;
        }
    }
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

auto get_interactions(t_type t, k_type k, v_type v) {
    auto cols = combinations(range(k), t);
    std::vector<std::vector<k_type>> col_sets;
    for (const auto& col_set : cols) {
        std::vector<k_type> to_add;
        for (const auto& new_col : col_set) {
            to_add.push_back(new_col);
        }
        col_sets.push_back(to_add);
    }
    std::vector<std::vector<v_type>> vals;
    if (t == 1) {
        for (v_type i = 0; i < v; i++) {
            std::vector<v_type> s{i};
            vals.push_back(s);
        }
    }
    else if (t == 2) {
        for (v_type i = 0; i < v; i++) {
            for (v_type j = 0; j < v; j++) {
                std::vector<v_type> s{i, j};
                vals.push_back(s);
            }
        }
    }
    else if (t == 3) {
        for (v_type i = 0; i < v; i++) {
            for (v_type j = 0; j < v; j++) {
                for (v_type k = 0; k < v; k++) {
                    std::vector<v_type> s{i, j, k};
                    vals.push_back(s);
                }
            }
        }
    }
    else if (t == 4) {
        for (v_type i = 0; i < v; i++) {
            for (v_type j = 0; j < v; j++) {
                for (v_type k = 0; k < v; k++) {
                    for (v_type l = 0; l < v; l++) {
                        std::vector<v_type> s{i, j, k, l};
                        vals.push_back(s);
                    }
                }
            }
        }
    }
    else if (t == 5) {
        for (v_type i = 0; i < v; i++) {
            for (v_type j = 0; j < v; j++) {
                for (v_type k = 0; k < v; k++) {
                    for (v_type l = 0; l < v; l++) {
                        for (v_type m = 0; m < v; m++) {
                            std::vector<v_type> s{i, j, k, l, m};
                            vals.push_back(s);
                        }
                    }
                }
            }
        }
    }
    else if (t == 6) {
        for (v_type i = 0; i < v; i++) {
            for (v_type j = 0; j < v; j++) {
                for (v_type k = 0; k < v; k++) {
                    for (v_type l = 0; l < v; l++) {
                        for (v_type m = 0; m < v; m++) {
                            for (v_type n = 0; n < v; n++) {
                                std::vector<v_type> s{i, j, k, l, m, n};
                                vals.push_back(s);
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

    std::vector<interaction_type> interactions;
    for (const auto& col_set : col_sets) {
        for (const auto& val_set : vals) {
            interaction_type I = std::make_pair(col_set, val_set);
            interactions.push_back(I);
        }
    }
    return interactions;
}

auto find_non_locating_sets(ca_type& A, t_type t, k_type k, v_type v, lambda_type lambda, unsigned int d) {
    auto interactions = get_interactions(t, k, v);

    // get all at most d;
    std::vector<d_set_type> d_sets;
    for (int i=1; i<=d; i++) {
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

auto read_ca_from_cagen(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    // skip first line as that does not have CA info
    std::getline(file, line);
    ca_type result;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<v_type> ca_line;
        for (int i; ss >> i;) {
            ca_line.push_back(i);
            if (ss.peek() == ',') {
                ss.ignore();
            }
        }
        result.push_back(ca_line);
    }
    return result;
}

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
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

int fitness(const ca_type& ind, d_type d, t_type t, k_type k, v_type v, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs) {
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
    }
    return score;
}

ca_type cross(const ca_type& p1, const ca_type& p2, d_type d, t_type t, k_type k, v_type v, lambda_type l) {
    int val = rand() % 2;
    int n = p1.size();
    ca_type child;
    if (val == 0) {
        auto rand_idx = rand() % p1.size();
        for (int i=0; i<rand_idx; i++) {
            child.push_back(p1[i]);
        }
        for (int i=rand_idx; i<n; i++) {
            child.push_back(p2[i]);
        }
    } else if (val == 1) {
        auto rand_idx1 = rand() % (p1.size());
        auto rand_idx2 = rand() % (p1.size());
        while (rand_idx1 == rand_idx2) {
            rand_idx1 = rand() % (p1.size());
            rand_idx2 = rand() % (p1.size());
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
    }
    return child;
}

ca_type mutate(const ca_type& p1, d_type d, t_type t, k_type k, v_type v, lambda_type l) {
    int val = rand() % 3;
    int n = p1.size();
    ca_type child = p1;
    if (val == 0) {
        auto rand_row = rand() % n;
        for (int col=0; col<p1[0].size(); col++) {
            child[rand_row][col] = rand() % v;
        }
    } else if (val == 1) {
        auto rand_col = rand() % k;
        for (int row=0; row<n; row++) {
            child[row][rand_col] = rand() % v;
        }
    } else {
        auto rand_row = rand() % p1.size();
        auto rand_col = rand() % p1[0].size();
        auto rand_val = rand() % v;
        child[rand_row][rand_col] = rand_val;
    }
    return child;
}

struct Ind_NonRecompute_Fitness {
    ca_type A;
    int fitness;
};

ca_type try_N(N_type N, d_type d, t_type t, k_type k, v_type v, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs) {
    ca_type s;
    int pop_size = 100;
    int num_gens = 100;
    std::vector<Ind_NonRecompute_Fitness> pop(pop_size);
    for (auto& elem : pop) {
        elem.A = random_array(N, k, v);
        elem.fitness = -1;
    }
    auto best_fitness = std::numeric_limits<int>::min();
    auto max_possible_fitness = non_locating_pairs.size();

    for (int gen=0; gen<num_gens; gen++) {
        std::vector<std::pair<int, Ind_NonRecompute_Fitness>> fitnesses;
        for (auto& I : pop) {
            int f = I.fitness;
            if (f == -1) {
                f = fitness(I.A, d, t, k, v, l, non_locating_pairs);
                I.fitness = f; 
            }
            
            if (f == max_possible_fitness) {
                return I.A;
            }
            fitnesses.push_back(std::make_pair(f,I));
        }
        std::sort(fitnesses.begin(), fitnesses.end(), [](const auto& first, const auto& second) {
            return first.first < second.first;
        });

        // get half the elems
        std::vector<Ind_NonRecompute_Fitness> new_vec;
        for (int i=pop_size/2; i < pop_size; i++) {
            new_vec.push_back(fitnesses[i].second);
        }
        pop = new_vec;
        new_vec.clear();

        // create a sub population with half elements
        while (new_vec.size() < pop_size / 2) {
            auto idx1 = rand() % pop.size();
            auto idx2 = rand() % pop.size();

            const auto& p1 = pop[idx1];
            const auto& p2 = pop[idx2];

            auto cross_percent = rand() % 10;
            auto mut_percent = rand() % 10;
            if (cross_percent == 0 && mut_percent < 3) {
                auto new_ind = cross(p1.A,p2.A,d,t,k,v,l);
                new_ind = mutate(new_ind,d,t,k,v,l);
                Ind_NonRecompute_Fitness true_new_ind;
                true_new_ind.A = new_ind;
                true_new_ind.fitness = -1;
                new_vec.push_back(true_new_ind);
            } else if (cross_percent == 0) {
                auto new_ind = cross(p1.A,p2.A,d,t,k,v,l);
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

ca_type go(d_type d, t_type t, k_type k, v_type v, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs) {
    int N = 3;
    
    bool succ_first = true;
    ca_type result;
    while (true) {
        result = try_N(N, d, t, k, v, l, non_locating_pairs);
        if (succ_first && result.size() > 0) {
            return result;
        }
        if (result.size() > 0) {
            break;
        }
        N *= 2;
        succ_first = false;
    }
    int N_hi = N;
    int N_lo = N / 2;
    while (N_lo < N_hi) {
        int N_mid = (N_lo + N_hi) / 2;
        auto result2 = try_N(N_mid, d, t, k, v, l, non_locating_pairs);
        if (result2.size() > 0) {
            N_hi = N_mid;
            result = result2;
        } else {
            N_lo = N_mid + 1;
        }
    }
    return result;
}

int main(int argc, char** argv) {
    const bool LLL_instead_of_file = false;
    k_type ks[] = {10, 15};
    for (t_type t=2; t <= 2; t++) {
        for (auto k : ks) {
            for (d_type d=1; d <= 1; d++) {  
                for (v_type v=2; v <= 2; v++) {
                    //assert(d < v);
                    for (lambda_type lambda=3; lambda <= 3; lambda++) {
                        auto filename = "./evaluation/" + std::to_string(v) + "^" + std::to_string(k) + "-t" + std::to_string(t) + "_l" + std::to_string(lambda);
                        auto ext = ".csv";
                        std::cout << "------------d=" << std::to_string(d) << " " << filename << "------------\n";
                        ca_type A;

                        auto start = high_resolution_clock::now();

                        if (LLL_instead_of_file) {
                            A = LLL_gen(t,k,v,lambda);
                            std::cout << "LLL done with " << A.size() << " rows.\n";
                        } else {
                            A = read_ca_from_cagen(filename + ext);
                            std::cout << "Read file with " << A.size() << " rows.\n";
                        }
                        auto non_locating_pairs = find_non_locating_sets(A, t, k, v, lambda, d);
                        std::cout << "There were " << non_locating_pairs.size() << " non-locating pairs\n";
                        auto stop = high_resolution_clock::now();

                        auto smallest_GA_rows = std::numeric_limits<N_type>::max();
                        auto ga_total_time = std::numeric_limits<int>::max();
                        for (int ga_run = 0; ga_run < 1; ga_run++) {
                            auto ga_start_time = high_resolution_clock::now();
                            auto ga_rows = go(d,t,k,v,lambda,non_locating_pairs);
                            auto ga_end_time = high_resolution_clock::now();
                            if (ga_rows.size() < smallest_GA_rows) {
                                smallest_GA_rows = ga_rows.size();
                                ga_total_time = duration_cast<milliseconds>(ga_end_time-ga_start_time).count();
                            }
                        }
                        

                        
                        auto stage1_diff = duration_cast<milliseconds>(stop-start).count();
                        auto total_time = stage1_diff + ga_total_time;

                        
                        std::cout << "Total N=" << A.size() + smallest_GA_rows << ", time=" << total_time << "ms.\n";
                        // std::ofstream out(write_filename);
                        // out << std::to_string(diff.count()) << "ms\n";
                        // // std::vector<std::tuple<d_set_type, d_set_type, int>> to_return;
                        // for (auto&& [dset_1, dset_2, num] : non_locating_pairs) {
                        //  auto s1 = d_set_to_str(dset_1);
                        //  auto s2 = d_set_to_str(dset_2);
                        //  auto to_write = s1 + s2 + std::to_string(num) + "\n";
                        //  // std::cout << to_write << "\n";
                        //  out << to_write;
                        // }
                        // out.close();

                        // // write created CA out to file if using LLL
                        // if (LLL_instead_of_file) {
                        //     std::ofstream A_out("./LLL_outputs/t" + std::to_string(t) + "_k" + std::to_string(k) + "_v" + std::to_string(v) + "_l" + std::to_string(lambda) + ".txt");
                        //     for (const auto& row : A) {
                        //         for (const auto& elem : row) {
                        //             A_out << elem << " ";
                        //         }
                        //         A_out << "\n";
                        //     }
                        //     A_out.close();
                        // }
                    }
                }
            }
        }
    }
}
