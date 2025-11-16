/* NEW: ryandougherty/locatingarraygenerator/LocatingArrayGenerator-degryse/phase2_ce/phase2_ce.h */
#pragma once
#include "../LocAG.h"

/**
 * @brief Runs the deterministic Conditional Expectation (CE) algorithm
 * adapted for Locating Arrays.
 *
 * This algorithm is based on Algorithm 4 from Dougherty et al. (2022) 
 *. It is adapted from a *covering* problem to a *locating* * problem.
 *
 * Instead of minimizing the expected number of *uncovered interactions*,
 * this algorithm minimizes the expected number of *undistinguished
 * d-set pairs*[cite: 331]. It builds rows deterministically by choosing
 * the factor level that minimizes this expectation, calculated over
 * a target number of rows, N.
 *
 * @param array A pointer to the LocatingArray object. This object
 * contains all problem parameters (k, vs, lambda, etc.) and
 * the list of 'undistinguished_pairs' from Phase 1.
 * The function will add the new, deterministically-generated
 * rows to 'array->array'.
 */
void run_phase_2_ce(LocatingArray *array);