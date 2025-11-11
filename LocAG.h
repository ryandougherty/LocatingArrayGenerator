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

// 2. Include phase1.h
#include "phase1/phase1.h"

// 3. Include phase2.h, which defines PercentGAFitnessInd
#include "phase2/phase2.h"

// 4. Define the LocatingArray struct expected by phase2_ce
struct LocatingArray {
    ca_type array;      // The actual array data (std::vector<std::vector<...>>)
    k_type k;           // Number of columns
    vs_type vs;         // Vector of levels (for varied levels)
    v_type v;           // Number of levels (assuming uniform)
    t_type t;           // Strength
    lambda_type lambda; // Lambda

    // This is the map of uncovered interactions that phase2_ce expects
    std::map<interaction_type, int> uncovered_interactions;
};


// 5. Finally, include phase2_ce.h, which *uses* LocatingArray
#include "phase2_ce/phase2_ce.h"


// --- Type Aliases for Chrono (from LocAG.cpp) ---
using high_resolution_clock = std::chrono::high_resolution_clock;
using milliseconds = std::chrono::milliseconds;

#endif // LOCAG_H