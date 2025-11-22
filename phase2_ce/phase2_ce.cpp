/* MODIFIED: ryandougherty/locatingarraygenerator/LocatingArrayGenerator-degryse/phase2_ce/phase2_ce.cpp */

#include "phase2_ce.h"
#include "../utils/utils.h" // For d_set_to_str
#include <cmath>
#include <limits>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <numeric> // For std::iota
#include <iomanip> // For std::setprecision

// --- Helper Functions for Binomial Probabilities ---

// Pre-calculated log-factorials for performance
std::vector<double> log_factorial_cache;
void precompute_log_factorials(int n) {
    // [WARNING FIX] Cast n to size_t for comparison
    if (log_factorial_cache.size() >= static_cast<size_t>(n + 1)) return;
    log_factorial_cache.resize(n + 1);
    log_factorial_cache[0] = 0.0;
    for (int i = 1; i <= n; ++i) {
        log_factorial_cache[i] = log_factorial_cache[i - 1] + std::log((double)i);
    }
}

// Calculates log(nCk) using the precomputed cache
double log_nCr(int n, int r) {
    if (r < 0 || r > n) return -std::numeric_limits<double>::infinity(); // Log(0)
    // [WARNING FIX] Cast n to size_t for comparison
    if (log_factorial_cache.size() <= static_cast<size_t>(n)) {
        precompute_log_factorials(n);
    }
    return log_factorial_cache[n] - log_factorial_cache[r] - log_factorial_cache[n - r];
}

/**
 * @brief Calculates P(X <= k) for X ~ Binomial(n, p)
 * [PRECISION FIX] Uses long double for the summation.
 */
double binomial_cdf(int k, int n, double p) {
    if (k < 0) return 0.0;
    if (p == 0.0) return 1.0; // P(X=0) is 1. P(X<=k) is 1 for k>=0.
    if (p == 1.0) return (k >= n) ? 1.0 : 0.0; // P(X=n) is 1.
    
    long double cdf = 0.0L; // <-- Use long double for accumulator
    
    double log_p = std::log(p);
    double log_1_p = std::log(1.0 - p);
    
    for (int i = 0; i <= k; ++i) {
        if (i > n) break;
        double log_pmf = log_nCr(n, i) + (double)i * log_p + (double)(n - i) * log_1_p;
        cdf += std::exp(log_pmf); // <-- Add double to long double
    }
    return (double)cdf; // Cast back to double
}

// --- Helper Functions for CE Probability ---

// -1 in a row means "indeterminate" or "*"
const v_type INDETERMINATE = -1; 

/**
 * @brief [HELPER] Checks if a given *full* row covers a d-set.
 * (Copied from phase2_greedy.cpp for deterministic check)
 */
/* MODIFIED HELPER in phase2_ce.cpp */
static bool row_covers_d_set(const std::vector<v_type>& row, const d_set_type& d_set) {
    if (d_set.empty()) return false;
    
    for (const auto& interaction : d_set) {
        bool covers_interaction = true;
        if (interaction.first.empty()) continue;

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
        // [FIX] If ANY interaction is covered, the d-set is covered.
        if (covers_interaction) {
            return true;
        }
    }
    // [FIX] Only return false if NONE were covered
    return false;
}


double prob_covers_intersection(const std::vector<v_type>& partial_row, const d_set_type& d_set, const vs_type& vs) {
    double prob = 1.0;
    std::map<k_type, v_type> col_constraints;

    for (const auto& interaction : d_set) {
        for (size_t i = 0; i < interaction.first.size(); ++i) {
            k_type col = interaction.first[i];
            v_type val = interaction.second[i];
            if (partial_row[col] != INDETERMINATE && partial_row[col] != val) {
                return 0.0; // Impossible to cover
            }
            if (col_constraints.count(col) && col_constraints[col] != val) {
                 return 0.0; // Conflicting interaction
            }
            col_constraints[col] = val;
        }
    }

    for (const auto& [col, val] : col_constraints) {
        if (partial_row[col] == INDETERMINATE) {
            prob *= (1.0 / (double)vs[col]);
        }
    }
    return prob;
}

/**
 * @brief Calculates the probability that a random completion of a
 * partial row will cover a single d-set.
 */
double prob_covers_dset(const std::vector<v_type>& partial_row, const d_set_type& d_set, const vs_type& vs) {
    if (d_set.empty()) return 0.0;

    double total_prob = 0.0;
    int n = d_set.size();
    
    // Iterate through all non-empty subsets (1 to 2^n - 1)
    // Using bit manipulation for subsets
    for (int i = 1; i < (1 << n); ++i) {
        d_set_type subset;
        int set_bits = 0;
        
        for (int j = 0; j < n; ++j) {
            if ((i >> j) & 1) {
                subset.push_back(d_set[j]);
                set_bits++;
            }
        }

        double p_intersection = prob_covers_intersection(partial_row, subset, vs);

        // Inclusion-Exclusion: Add if odd size, Subtract if even size
        if (set_bits % 2 == 1) {
            total_prob += p_intersection;
        } else {
            total_prob -= p_intersection;
        }
    }
    
    // Clamp result for floating point errors
    if (total_prob < 0.0) return 0.0;
    if (total_prob > 1.0) return 1.0;
    return total_prob;
}

/**
 * @brief [COMPILER FIX] Calculates the probability that a random completion of
 * 'partial_row' will *distinguish* d_set1 and d_set2.
 * P(distinguish) = P(covers1) + P(covers2) - 2 * P(covers1 and covers2)
 */
/* MODIFIED in phase2_ce.cpp */
double prob_distinguishes(const std::vector<v_type>& partial_row, const d_set_type& d_set1, const d_set_type& d_set2, const vs_type& vs) {
    
    double p1 = prob_covers_dset(partial_row, d_set1, vs); // P(C1)
    double p2 = prob_covers_dset(partial_row, d_set2, vs); // P(C2)

    // Union of the two d-sets
    d_set_type d_union = d_set1;
    d_union.insert(d_union.end(), d_set2.begin(), d_set2.end());
    
    // With new logic, this calculates P(C1 OR C2)
    double p_union_or = prob_covers_dset(partial_row, d_union, vs);

    // Formula for symmetric difference using Union probability:
    // P(Diff) = 2 * P(Union) - P(A) - P(B)
    double result = 2.0 * p_union_or - p1 - p2;
    
    if (result < 0.0) return 0.0;
    if (result > 1.0) return 1.0;
    return result;
}


/**
 * @brief Calculates P(pair remains NOT lambda-distinguished)
 */
double prob_remains_undistinguished(int N, int M, double p, double cr, int lambda_remaining) {
    if (lambda_remaining <= 0) return 0.0; 
    if (p == 0.0 && cr == 0.0) return 1.0; 

    int N_remaining = N - M - 1;
    if (N_remaining < 0) N_remaining = 0;

    double p_fail_if_dist = binomial_cdf(lambda_remaining - 2, N_remaining, p);
    double p_fail_if_not_dist = binomial_cdf(lambda_remaining - 1, N_remaining, p);

    return (p_fail_if_dist * cr) + (p_fail_if_not_dist * (1.0 - cr));
}


/**
 * @brief Calculates the smallest N.
 * [PRECISION FIX] Uses long double for the summation.
 */
int calculate_N(const LocatingArray* array, const std::map<size_t, double>& p_dist_cache, int M, int N_start, const std::set<size_t>& remaining_pair_indices) {
    int N = N_start;
    precompute_log_factorials(N + 100); 

    while (true) {
        long double total_expected_failures = 0.0L; // <-- Use long double
        
        for (const auto& index : remaining_pair_indices) {
            const auto& [d_set1, d_set2, times_separated] = array->undistinguished_pairs[index];
            int lambda_remaining = array->lambda - times_separated;

            double p = p_dist_cache.at(index);
            total_expected_failures += binomial_cdf(lambda_remaining - 1, N - M, p);
        }

        if (total_expected_failures < 1.0L) { // <-- Compare to long double
            break; // Found our N
        }
        N++; // Try a larger N
        if (N > N_start + 2000) { 
             std::cout << "  (CE) WARNING: Could not find N < " << N << ". Using " << N << "." << std::endl;
             break;
        }
        // [WARNING FIX] Cast N to size_t for comparison
        if (static_cast<size_t>(N) > log_factorial_cache.size() - 10) {
            precompute_log_factorials(N + 100);
        }
    }
    return N;
}

/**
 * @brief Main CE algorithm implementation.
 * [PRECISION FIX] Uses long double for the summation.
 */
void run_phase_2_ce(LocatingArray *array) {

    int k = array->k;
    const vs_type& vs = array->vs;
    int M = array->array.size(); 
    int lambda = array->lambda;

    std::set<size_t> remaining_pair_indices;
    for(size_t i = 0; i < array->undistinguished_pairs.size(); ++i) {
        if (std::get<2>(array->undistinguished_pairs[i]) < lambda) {
            remaining_pair_indices.insert(i);
        }
    }

    std::map<size_t, double> p_dist_cache;
    std::vector<v_type> empty_row(k, INDETERMINATE);
    for (size_t i = 0; i < array->undistinguished_pairs.size(); ++i) {
        const auto& [d_set1, d_set2, ts] = array->undistinguished_pairs[i];
        p_dist_cache[i] = prob_distinguishes(empty_row, d_set1, d_set2, vs);
    }
    
    int N_target = calculate_N(array, p_dist_cache, M, M + 1, remaining_pair_indices);
    std::cout << "  (CE) Initial N_target calculated as: " << N_target << " (M=" << M << ")" << std::endl;

    while (!remaining_pair_indices.empty()) {
        std::vector<v_type> new_row(k, INDETERMINATE);
        
        for (int j = 0; j < k; ++j) { 
            long double min_expected_failures = std::numeric_limits<long double>::infinity(); // <-- Use long double
            v_type best_v = 0;

            for (v_type v = 0; v < vs[j]; ++v) { 
                new_row[j] = v; 
                long double current_expected_failures = 0.0L; // <-- Use long double

                for (const auto& index : remaining_pair_indices) {
                    const auto& [d1, d2, ts] = array->undistinguished_pairs[index];
                    int lambda_remaining = lambda - ts;
                    
                    double cr = prob_distinguishes(new_row, d1, d2, vs);
                    double p = p_dist_cache.at(index);
                    
                    current_expected_failures += prob_remains_undistinguished(
                        N_target, M, p, cr, lambda_remaining);
                }

                if (current_expected_failures < min_expected_failures) {
                    min_expected_failures = current_expected_failures;
                    best_v = v;
                }
            } 
            
            new_row[j] = best_v;
        } 

        array->array.push_back(new_row);
        M++; 

        std::vector<size_t> distinguished_this_round;
        for (const auto& index : remaining_pair_indices) {
            auto& pair_tuple = array->undistinguished_pairs[index];
            auto& times_separated = std::get<2>(pair_tuple);
            
            const auto& d_set1 = std::get<0>(pair_tuple);
            const auto& d_set2 = std::get<1>(pair_tuple);

            // [LOGIC FIX] Use the deterministic checker
            bool covers1 = row_covers_d_set(new_row, d_set1);
            bool covers2 = row_covers_d_set(new_row, d_set2);

            if ((covers1 && !covers2) || (!covers1 && covers2)) {
                times_separated++;
                if (times_separated >= lambda) {
                    distinguished_this_round.push_back(index);
                }
            }
        }

        for (const auto& index : distinguished_this_round) {
            remaining_pair_indices.erase(index);
        }
        
        std::cout << "  (CE) Row " << M << " built. "
                  << remaining_pair_indices.size() << " pairs remaining to distinguish." << std::endl;

        if (!remaining_pair_indices.empty()) {
            int old_N_target = N_target;
            N_target = calculate_N(array, p_dist_cache, M, N_target, remaining_pair_indices); 
            if (N_target != old_N_target) {
                std::cout << "  (CE) N_target recalculated as: " << N_target << std::endl;
            }
        }

        if (M > N_target + k*20 && M > 500) { 
             std::cout << "  (CE) WARNING: Algorithm seems stuck (M > N_target). Forcefully exiting." << std::endl;
             std::cout << "  (CE) " << remaining_pair_indices.size() << " pairs were left undistinguished." << std::endl;
             break;
        }
    } 
}