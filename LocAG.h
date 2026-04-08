/* MODIFIED: ryandougherty/locatingarraygenerator/LocatingArrayGenerator-degryse/LocAG.h */
#ifndef LOCAG_H
#define LOCAG_H

// --- Standard Libraries (from LocAG.cpp) ---
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

// --- Third-party (from LocAG.cpp) ---
#include "combinations.hpp"
#include "enumerate.hpp"
#include "product.hpp"
#include "range.hpp"
#include "robin_hood.h"
#include "zip.hpp"

// --- Project Headers (in correct order) ---

// 1. Include utils.h first, as it defines the basic types
// (ca_type, t_type, interaction_type, etc.) and hashers.
#include "utils/utils.h"

// 1.5 Include the interaction codec for integer encoding of interactions/d-sets.
#include "utils/interaction_codec.h"

// 2. Include phase1.h
#include "phase1/phase1.h"

// 3. Include phase2.h, which defines PercentGAFitnessInd
#include "phase2/phase2.h"

// 4. Define the LocatingArray struct
// This struct holds the data for both the greedy and CE locating problem.
struct LocatingArray {
    ca_type array;      // The full array (initial + new rows)
    k_type k;           // Number of columns
    vs_type vs;         // Vector of levels (for varied levels)
    v_type v;           // Average number of levels (for heuristic)
    t_type t;           // Strength
    lambda_type lambda; // Lambda
    d_type d;           // d-set size
    bool is_detecting;  // Flag for detecting vs locating

    // Codec for encoding/decoding interactions and d-sets as integers.
    InteractionCodec codec;

    // The list of pairs that Phase 1 found to be undistinguished.
    // Each entry is (d_set_id_1, d_set_id_2, times_separated).
    // Use codec.decode_d_set(id) to recover the interaction_type objects.
    std::vector<undist_pair_type> undistinguished_pairs;
};


// 5. Finally, include the headers for the phase 2 algorithms
//    (Renamed old ce -> greedy)
#include "phase2_greedy/phase2_greedy.h"
//    (This is the new, true CE algorithm)
#include "phase2_ce/phase2_ce.h"


// --- Type Aliases for Chrono (from LocAG.cpp) ---
using high_resolution_clock = std::chrono::high_resolution_clock;
using milliseconds = std::chrono::milliseconds;

#endif // LOCAG_H