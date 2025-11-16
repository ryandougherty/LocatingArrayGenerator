/* ----------------------------------------------------------------------------
 * LocAG.cpp
 *
 * This is the main entry point for the Locating Array Generator (LocAG) tool.
 *
 * The program operates in two main stages:
 * 1.  **Phase 1 (Analysis):** Reads an existing covering array (CA) from a file,
 * identifies all "non-locating" or "non-detecting" d-sets. These are
 * sets of interactions that are not sufficiently distinguished by the
 * rows of the array.
 * 2.  **Phase 2 (Genetic Algorithm):** Uses a Genetic Algorithm (GA) to
 * evolve a set of *new* rows. These new rows are added to the
 * original array to "fix" the non-valid pairs found in Phase 1,
 * resulting in a valid locating or detecting array.
 *
 * Usage: ./LocAG <config_name> <array_type> <execution_policy>
 * - config_name: Name of a pre-defined configuration (e.g., "Apache") or
 * a parameter string (e.g., "2^30").
 * - array_type: "locating" or "detecting".
 * - execution_policy: "serial" or "parallel" (for the GA).
 * ----------------------------------------------------------------------------
 */

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

// --- Third-party libraries ---
#include "combinations.hpp" // For iterating through combinations (e.g., columns)
#include "enumerate.hpp"   // For iterating with (index, value) pairs
//#include "flat_hash_map.hpp" // Alternate hash map
#include "product.hpp"     // For Cartesian products
#include "range.hpp"       // For Python-style range()
#include "robin_hood.h"    // High-performance hash map
#include "zip.hpp"         // For iterating over multiple containers in parallel

// --- Project-specific headers ---
#include "LocAG.h"         // Type definitions (ca_type, t_type, etc.)
#include "utils/utils.h"   // Utility functions (printing, random numbers)
#include "phase1/phase1.h" // Phase 1: Finding non-valid pairs
#include "phase2/phase2.h" // Phase 2: Genetic algorithm for fixing pairs
#include "phase2_greedy/phase2_greedy.h"   // For the new CE algorithm
#include "phase2_ce/phase2_ce.h"   // For the new CE algorithm

/**
 * @brief Checks a covering array 'A' for t-way coverage with a given lambda.
 *
 * This function iterates through all t-combinations of columns. For each
 * combination, it counts the occurrences of every possible value-tuple.
 * If any tuple appears less than 'lambda' times, or if not all
 * v^t tuples are present, it returns the columns that failed the check.
 *
 * @param A The covering array (ca_type).
 * @param t The strength (t-way coverage).
 * @param k The number of columns.
 * @param vs The vector of levels for each column.
 * @param lambda The minimum number of times each t-tuple must appear.
 * @return A vector of column indices that failed coverage, or an empty
 * vector if 'A' is a valid (t, k, v, lambda) CA.
 */
auto first_uncovered_cols(ca_type A, const t_type t, const k_type k, const vs_type& vs, const lambda_type lambda) {
    std::vector<v_type> row_in_A(t, 0);
    std::vector<k_type> cols_to_return;
    
    // Iterate over every combination of 't' columns
    for (const auto& cols : combinations(range(k), t)) {
        // Use a hash map to count occurrences of each value-tuple
        robin_hood::unordered_flat_map<std::vector<v_type>, int, VectorHasher> c;
        
        // Project each row onto the selected columns
        for (const auto& row : A) {
            for (size_t i = 0; i < cols.size(); i++) {
                row_in_A[i] = row[cols[i]];
            }
            // Increment the count for this specific tuple
            if (c.find(row_in_A) != c.end()) {
                c[row_in_A] += 1;
            }
            else {
                c[row_in_A] = 1;
            }
        }

        // Check if all possible interactions are present
        // --- FIX: This should be based on the product of levels for the chosen columns ---
        long long expected_interactions = 1;
        for(k_type col_idx : cols) {
            expected_interactions *= vs[col_idx];
        }

        if (c.size() != (size_t)expected_interactions) {
        // --- END FIX ---
            for (const auto& col : cols) {
                cols_to_return.push_back(col);
            }
            return cols_to_return; // Return the failing columns
        }
        else {
            // If all interactions are present, check if they meet the lambda requirement
            const auto& it = *std::min_element(std::begin(c), std::end(c),
                [](const auto& l, const auto& r) { return l.second < r.second; });
            if (it.second < lambda) {
                for (const auto& col : cols) {
                    cols_to_return.push_back(col);
                }
                return cols_to_return; // Return the failing columns
            }
        }
    }
    return std::vector<k_type>(); // All columns are covered
}

/**
 * @brief Reads a covering array from a CSV file (cagen format).
 *
 * This parser handles a specific format where:
 * - The first line is skipped (assumed to be a header).
 * - Each subsequent line is a row of the array.
 * - Values are comma-separated.
 * - A "*" (asterisk) is treated as a "don't care" or wildcard, and a
 * randomly chosen valid value is substituted.
 *
 * @param filename The path to the .csv file.
 * @param vs The vector of levels for each column (used for wildcard generation).
 * @return The parsed covering array (ca_type).
 */
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
                // If wildcard, substitute a random value from 0 to vs[col_idx]-1
                ca_line.push_back(any_int(rng) % vs[col_idx]);
            } else {
                // Otherwise, parse the integer value
                ca_line.push_back(std::stoi(elem));
            }
        }
        result.push_back(ca_line);
    }
    return result;
}

/**
 * @brief Memoization helper to get the rows for an interaction.
 *
 * Checks if the rows for a given 'interaction' are already in 'rows_map'.
 * If yes, returns the cached set of rows.
 * If no, computes the rows using 'rows_of_interaction', stores them in
 * the map, and then returns them.
 *
 * @param rows_map A map caching interactions and their corresponding row sets.
 * @param interaction The interaction to look up.
 * @param ind The array to search within.
 * @return A (const) set of row indices where the interaction appears.
 */
auto lookup_or_assign_interaction_map(std::map<interaction_type, robin_hood::unordered_flat_set<int>>& rows_map, const interaction_type& interaction, const ca_type& ind) {
    if (rows_map.find(interaction) != rows_map.end()) {
        // Found in cache
        const auto& rows = rows_map[interaction];
        return rows;
    } else {
        // Not in cache, compute and store
        const auto& rows = rows_of_interaction(interaction, ind);
        rows_map[interaction] = rows;
        return rows;
    } 
}

/**
 * @brief Parses a vector of "v^k" strings into a single 'vs_type' vector.
 *
 * Example: {"2^3", "4^1"} becomes {2, 2, 2, 4}
 *
 * @param exp_params A vector of strings (e.g., "2^158").
 * @return A 'vs_type' (vector<v_type>) listing the levels for all columns.
 */
auto parse_vs(const std::vector<std::string>& exp_params) {
    vs_type new_params;
    for (const auto& str : exp_params) {
        auto idx = str.find("^");
        // Parse 'v' (the number of levels)
        v_type v = std::atoi(std::string(str.begin(), str.begin() + idx).c_str());
        // Parse 'k' (the number of columns with this many levels)
        auto k = std::atoi(std::string(str.begin() + idx+1, str.end()).c_str());
        // Add 'v' to the result vector 'k' times
        for (int i=0; i<k; i++) {
            new_params.push_back(v);
        }
    }
    return new_params;
}

// A hardcoded map of well-known software/system configurations
// used for generating locating/detecting arrays.
// The strings define the parameters (levels^columns)
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
    //other papers with these for comparison***
};

/**
 * @brief Looks up a configuration by name and generates parameters.
 *
 * If 'config_name' is found in the global 'configs' map, it parses
 * those parameters.
 * If not found, it assumes 'config_name' *is* the parameter string
 * itself (e.g., "2^30").
 *
 * @param config_name The name to look up (e.g., "Apache") or a param string.
 * @param t The strength parameter.
 * @param lambda The lambda parameter.
 * @return A pair containing:
 * 1. The 'vs_type' vector of column levels.
 * 2. The generated filename for the input array.
 */
auto lookup_config_and_params(const std::string& config_name, const t_type t, const lambda_type lambda) {

    std::string prefix = "./evaluation/";
    vs_type the_vals;
    if (configs.find(config_name) != configs.end()) {
        // It's a known config name
        the_vals = parse_vs(configs.at(config_name));
        auto filename = prefix + config_name + "/" + config_name + "_" + std::to_string(t) + "_" + std::to_string(lambda) + ".csv";
        return std::make_pair(the_vals, filename);
    } else {
        // Not a known name, assume it's a parameter string like "2^30"
        std::vector<std::string> the_new_config{config_name};
        the_vals = parse_vs(the_new_config);
        auto filename = prefix + config_name + "-t" + std::to_string(t) + "_l" + std::to_string(lambda) + ".csv";
        return std::make_pair(the_vals, filename);
    }
    
}

/**
 * @brief Main function for the Locating Array Generator.
 */
int main(int argc, char** argv) {

    // --- Argument Parsing ---
    if (argc != 5) {
        std::cerr << "Usage: ./LocAG <name of config> <array_type> <method> <execution_policy>\n";
        std::cerr << "array_type can be 'locating' or 'detecting'\n";
        std::cerr << "method can be 'ga', 'greedy', or 'ce'\n";
        std::cerr << "execution_policy can be 'serial' or 'parallel'\n";
        return -1;
    }

    std::string algorithm_type = "ga"; // Default
    std::string array_type = argv[2];
    const std::string policy = argv[4]; // "serial" or "parallel"
    
    // --- MODIFIED: Updated algorithm selection ---
    if (strcmp(argv[3],"greedy") == 0) {
        algorithm_type = "greedy";
    }
    else if (strcmp(argv[3], "ce") == 0) {
        algorithm_type = "ce";
    }
    else if (strcmp(argv[3], "ga") != 0) {
        std::cerr << "Invalid method: " << argv[3] << ". Must be 'ga', 'greedy', or 'ce'." << std::endl;
        return -1;
    }
    

    // --- Main Experiment Loop ---
    // Note: These loops are hardcoded to only run d=1, t=2
    for (d_type d = 1; d <= 1; d++) {
        for (t_type t = 2; t <= 2; t++) {
            for (lambda_type lambda = 1; lambda <= 4; lambda++) {
                const std::string config_name = argv[1];
                
                // Get parameters (vs) and input filename
                const auto& [vs, filename] = lookup_config_and_params(config_name, t, lambda);

                // --- Set Experiment Flags ---
                const bool d_bar = true; // Use "at most d" interactions
                const bool t_bar = true; // Use "at most t" interactions
                bool is_detecting = false;
                const int X = 10; // Number of partitions for Phase 1 heuristic

                // 'd' must be less than the smallest column level
                assert(d < *std::min_element(vs.begin(), vs.end()));

                std::cout << "------------d=" << std::to_string(d) << ", t=" << std::to_string(t) << ", lambda=" << std::to_string(lambda) << ", filename=" << filename << "------------\n";

                // --- STAGE 1: Analysis ---
                auto start = high_resolution_clock::now();
                
                // Read the initial covering array
                ca_type A = read_ca_from_cagen(filename, vs);
                // std::cout << "Read file with " << A.size() << " rows.\n";

                // Find all pairs of d-sets that are not correctly located/detected
                std::vector<std::tuple<d_set_type, d_set_type, int>> non_valid_pairs;
                if (array_type == "locating") {
                    non_valid_pairs = find_non_locating_sets(A, t, vs, lambda, d, d_bar, t_bar, X);
                }
                else if (array_type == "detecting") {
                    is_detecting = true;
                    non_valid_pairs = find_non_detecting_sets(A, t, vs, lambda, d, d_bar, t_bar, X);
                }
                else {
                    std::cerr << "Array type " + array_type + " is not valid.\n";
                    return -1;
                }
                auto stop = high_resolution_clock::now();
                auto first_stage_N = A.size();
                auto first_stage_time = duration_cast<milliseconds>(stop-start).count();

                std::cout << "First Stage N=" << first_stage_N << ", Time=" << first_stage_time << "\n"; 
                std::cout << "There are " << non_valid_pairs.size() << " remaining non-locating pairs\n";

                /* ------------------------------------- Stage 2: GA / Density ------------------------------ */

                // If the initial array is already valid, we're done.
                if (non_valid_pairs.size() == 0) {
                    // Create a filename for the output, even if no new rows were added
                    std::string output_filename = "./" + config_name + "_t" + std::to_string(t) + "_l" + std::to_string(lambda) + "_final.csv";
                    save_array_to_csv(A, output_filename); // Save the original array
                    std::cout << "Array is already valid. Saved initial array to " << output_filename << "\n";
                    break;
                }

                std::vector<PercentGAFitnessInd> pareto;

                if (algorithm_type == "ga") {
                    std::cout << "--- Running GA Algorithm (phase2) ---" << std::endl;
                    // Run the genetic algorithm (Phase 2) to find new rows
                    // This GA optimizes the *percentage* of pairs to fix at each step
                    pareto = percent_GA(d,t,vs,lambda,non_valid_pairs,false, is_detecting, policy);
                    
                }
                // --- GLUE CODE FOR DENSITY ALGORITHM ---
                else if (algorithm_type == "greedy") {
                    std::cout << "--- Running Density Algorithm (phase2_ce) for Locating ---" << std::endl;

                    // 1. Create and populate the LocatingArray struct
                    LocatingArray loc_array;
                    loc_array.array = A; // The initial array from Phase 1
                    loc_array.vs = vs;
                    loc_array.t = t;
                    loc_array.lambda = lambda;
                    loc_array.k = vs.size();
                    loc_array.d = d;
                    loc_array.is_detecting = is_detecting;
                    
                    if (!vs.empty()) {
                        double sum_levels = 0.0;
                        for(v_type level : vs) { sum_levels += level; }
                        loc_array.v = static_cast<v_type>(std::round(sum_levels / vs.size())); 
                    } else {
                        loc_array.v = 0; // Should not happen given configs
                    }

                    // 2. Assign the non-valid pairs to the struct.
                    // The 'density' algorithm will now operate on this list.
                    loc_array.undistinguished_pairs = non_valid_pairs;
                    
                    // 3. Start timer and call the function
                    auto ce_start = high_resolution_clock::now();
                    run_phase_2_greedy(&loc_array); // Pass a pointer
                    auto ce_stop = high_resolution_clock::now();
                    auto ce_time = duration_cast<milliseconds>(ce_stop - ce_start).count();

                    // 4. Extract new rows. The function modified loc_array.array directly.
                    ca_type new_rows;
                    if (loc_array.array.size() > first_stage_N) {
                        new_rows.insert(new_rows.end(),
                                        loc_array.array.begin() + first_stage_N,
                                        loc_array.array.end());
                    }

                    // 5. Manually create a single result and add it to 'pareto'
                    PercentGAFitnessInd density_result;
                    density_result.N = new_rows.size();
                    density_result.time = ce_time;
                    density_result.generated_rows = new_rows;
                    density_result.percents = {1.0}; // 'percents' isn't applicable, so use a placeholder

                    pareto.push_back(density_result);
                }
                // --- END GLUE CODE ---
                // --- NEW: Added block for "ce" ---
                else if (algorithm_type == "ce") {
                    std::cout << "--- Running Conditional Expectation Algorithm (phase2_ce) for Locating ---" << std::endl;

                    // 1. Create and populate the LocatingArray struct
                    LocatingArray loc_array;
                    loc_array.array = A; // The initial array from Phase 1
                    loc_array.vs = vs;
                    loc_array.t = t;
                    loc_array.lambda = lambda;
                    loc_array.k = vs.size();
                    loc_array.d = d;
                    loc_array.is_detecting = is_detecting;
                    
                    if (!vs.empty()) {
                        double sum_levels = 0.0;
                        for(v_type level : vs) { sum_levels += level; }
                        loc_array.v = static_cast<v_type>(std::round(sum_levels / vs.size())); 
                    } else {
                        loc_array.v = 0; // Should not happen given configs
                    }

                    // 2. Assign the non-valid pairs to the struct.
                    loc_array.undistinguished_pairs = non_valid_pairs;
                    
                    // 3. Start timer and call the new CE function
                    auto ce_start = high_resolution_clock::now();
                    run_phase_2_ce(&loc_array); // <-- NEW function call
                    auto ce_stop = high_resolution_clock::now();
                    auto ce_time = duration_cast<milliseconds>(ce_stop - ce_start).count();

                    // 4. Extract new rows.
                    ca_type new_rows;
                    if (loc_array.array.size() > first_stage_N) {
                        new_rows.insert(new_rows.end(),
                                        loc_array.array.begin() + first_stage_N,
                                        loc_array.array.end());
                    }

                    // 5. Manually create a single result and add it to 'pareto'
                    PercentGAFitnessInd ce_result;
                    ce_result.N = new_rows.size();
                    ce_result.time = ce_time;
                    ce_result.generated_rows = new_rows;
                    ce_result.percents = {1.0}; 

                    pareto.push_back(ce_result);
                }
                // --- END NEW BLOCK ---

                // Print the results from the Pareto front
                // Each result is a trade-off between (total rows) and (total time)
                // --- Save the results from the Pareto front ---
                int pareto_solution_index = 0;
                std::cout << "--- Pareto Front Solutions ---" << std::endl;
                for (const auto& ind : pareto) {
                    // 1. Combine the initial array 'A' with the new rows
                    ca_type final_array = A;
                    final_array.insert(final_array.end(), ind.generated_rows.begin(), ind.generated_rows.end());

                    // 2. Create a unique output filename
                    std::string output_filename = "./results/" + config_name + "_t" + std::to_string(t) + "_l" + 
                                                std::to_string(lambda) + "_pareto_" + 
                                                std::to_string(pareto_solution_index) + ".csv";
                    
                    // 3. Save the final array
                    save_array_to_csv(final_array, output_filename);

                    // 4. Print results to console
                    std::cout << "Solution " << pareto_solution_index << ": N total=" << final_array.size() 
                            << " (N_initial=" << first_stage_N << ", N_added=" << ind.N << ")"
                            << ", Time total=" << first_stage_time + ind.time << " ms"
                            << ", Saved to: " << output_filename << ", percents=";
                    print_vec(ind.percents);
                    std::cout << "\n";
                    
                    pareto_solution_index++;
                }
            } // End lambda loop
        } // End t loop
    } // End d loop
    
    return 0;
} // End main

// --- Research/Todo Comments ---
//Do checks in parallel and for main loop
//Do GA double for loop in parallel
//Experimenting (look at papers on CA and software testing) look in transactions in sofwater engineering for research questions
//Research Questions: 
//1. Is this tool faster than previous tools? Is it more memory efficient? Does it run better in parallel? 
//Start writing a summary of the method (based on other paper) and new work