/* ----------------------------------------------------------------------------
 * phase2.cpp
 *
 * Implements Phase 2 of the LocAG algorithm: the optimization stage.
 *
 * This file contains two main optimization algorithms:
 *
 * 1.  A "meta" Genetic Algorithm (GA) (`percent_GA`):
 * - The "individuals" in this GA are not arrays, but *strategies*.
 * - A strategy is a vector of percentages (e.g., {0.02, 0.1, 0.5, 1.0}).
 * - This strategy means "First, generate rows to fix 2% of pairs, then
 * generate new rows to fix 10% of remaining, etc., until 100% are fixed."
 * - The GA's goal is to find a strategy (a `percents` vector) that
 * minimizes both the *total number of rows* (N) and the *total time*.
 * - This GA can be run in parallel, where each thread evaluates a
 * different strategy.
 *
 * 2.  A Simulated Annealing (SA) algorithm (`try_N_SA`):
 * - This is the "inner" search algorithm used by the GA to execute
 * one step of its strategy (e.g., "fix 2% of pairs").
 * - Given a target number of rows (N) and a list of pairs to fix,
 * the SA tries to find a *single* array of size N that fixes them all.
 * - It uses a binary search (`go` function) to find the *minimal* N
 * that the SA can successfully find a solution for.
 * ----------------------------------------------------------------------------
 */

#include "../utils/utils.h"
#include "phase2.h"
#include "../phase1/phase1.h" // For rows_of_interaction
#include <execution> // For std::execution::par

/**
 * @brief Fitness function: Calculates how many pairs an array 'ind' separates.
 *
 * This function is the core of the Simulated Annealing. It checks an
 * individual 'ind' (a candidate array) against a list of 'non_locating_pairs'.
 *
 * @param ind The candidate array (an "individual") to evaluate.
 * @param non_locating_pairs The list of (D1, D2, count) pairs to fix.
 * @param threshold An optimization: if the score reaches this, stop early.
 * @param is_detecting Flag to change requirement from lambda to 1.
 * @return The score (number of pairs successfully separated).
 */
int fitness(const ca_type& ind, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, const int& threshold, bool is_detecting) {
    int score = 0;
    
    // Memoization map to cache row sets for d-sets *within this individual*
    std::unordered_map<d_set_type, std::vector<N_type>, DSetHasher> rows_map;
    
    // Lambda helper to get rows for a d-set, using the cache
    auto rows_of_dset = [=,&rows_map](const d_set_type& d_set) {
        if (rows_map.find(d_set) != rows_map.end()) {
            return rows_map[d_set];
        } else {
            // Not in cache, compute it
            robin_hood::unordered_set<N_type> the_rows;
            for (const auto& interaction : d_set) {
                // Find rows *in the new array 'ind'*
                const auto& rows = rows_of_interaction(interaction,ind);
                the_rows.insert(rows.begin(), rows.end());
            }
            // Store as a sorted vector for faster set operations
            std::vector<int> vrows(the_rows.begin(), the_rows.end());
            std::sort(vrows.begin(), vrows.end());
            rows_map[d_set] = vrows;
            return vrows;
        }
    };

    // Check every pair that needs fixing
    for (const auto& [dset_1, dset_2, num_times_sep_already] : non_locating_pairs) {
        // Calculate how many *more* separating rows are needed
        auto requirement = (is_detecting ? 1 : l) - num_times_sep_already;
        
        // Get row sets from the new array 'ind'
        auto rows1 = rows_of_dset(dset_1);
        auto rows2 = rows_of_dset(dset_2);
        
        // Calculate the symmetric difference (number of separating rows)
        int n = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());
        
        if (n >= requirement) {
            score += 1; // This pair is now fixed
        }

        if (score >= threshold) {
            return threshold + 1; // Early exit optimization
        }
    }
    return score;
}

/**
 * @brief GA Crossover operator (for the deprecated 'try_N' GA).
 *
 * Performs one-point or two-point crossover on the *rows* of two parent arrays.
 *
 * @param p1 Parent 1 array.
 * @param p2 Parent 2 array.
 * @return A new child array.
 */
ca_type cross(const ca_type& p1, const ca_type& p2, d_type d, t_type t, const vs_type& vs, lambda_type l, std::mt19937& rng) {
    int val = any_int(rng) % 2;
    int n = p1.size();
    ca_type child;
    if (val == 0) {
        // One-point crossover
        auto rand_idx = any_int(rng) % p1.size();
        for (int i=0; i<rand_idx; i++) {
            child.push_back(p1[i]);
        }
        for (int i=rand_idx; i<n; i++) {
            child.push_back(p2[i]);
        }
    } else if (p1.size() != 1) {
        // Two-point crossover
        auto rand_idx1 = any_int(rng) % (p1.size());
        auto rand_idx2 = any_int(rng) % (p1.size());
        while (rand_idx1 == rand_idx2) {
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
        // Failsafe for N=1
        val = any_int(rng) % 2;
        if (val == 0) {
            child = p1;
        } else {
            child = p2;
        }
    }
    return child;
}

/**
 * @brief Mutation operator (used by Simulated Annealing).
 *
 * Randomly applies one of three mutations:
 * 1. Mutate a single cell.
 * 2. Mutate an entire row.
 * 3. Mutate an entire column.
 *
 * @param p1 The array to mutate.
 * @return A new, mutated array.
 */
ca_type mutate(const ca_type& p1, d_type d, t_type t, const vs_type& vs, lambda_type l, std::mt19937& rng) {
    int val = any_int(rng) % 3;
    int n = p1.size();
    ca_type child = p1;
    if (val == 0) {
        // Mutate an entire row
        auto rand_row = any_int(rng) % n;
        for (int col=0; col<p1[0].size(); col++) {
            child[rand_row][col] = any_int(rng) % vs[col];
        }
    } else if (val == 1) {
        // Mutate an entire column
        auto rand_col = any_int(rng) % vs.size();
        for (int row=0; row<n; row++) {
            child[row][rand_col] = any_int(rng) % vs[rand_col];
        }
    } else {
        // Mutate a single cell
        auto rand_row = any_int(rng) % p1.size();
        auto rand_col = any_int(rng) % p1[0].size();
        auto rand_val = any_int(rng) % vs[rand_col];
        child[rand_row][rand_col] = rand_val;
    }
    return child;
}

// Struct to hold a GA individual (array) and its cached fitness
struct Ind_NonRecompute_Fitness {
    ca_type A;
    int fitness;
};

/**
 * @brief A standard Genetic Algorithm (deprecated, not used by main).
 *
 * Tries to find an array of size 'N' that satisfies a 'percent' of
 * the 'non_locating_pairs'.
 *
 * @return The successful array, or an empty array on failure.
 */
ca_type try_N(N_type N, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, double percent, bool is_detecting, std::mt19937& rng) {

    ca_type s;
    int pop_size = 100;
    int num_gens = 50;

    // Initialize population
    std::vector<Ind_NonRecompute_Fitness> pop(pop_size);
    for (auto& elem : pop) {
        elem.A = random_array(N, vs.size(), vs);
        elem.fitness = -1; // -1 means fitness not yet computed
    }

    const int max_possible_fitness = non_locating_pairs.size() * percent;

    for (int gen=0; gen<num_gens; gen++) {
        // Calculate fitness for all individuals
        std::vector<std::pair<int, Ind_NonRecompute_Fitness>> fitnesses;
        for (auto& I : pop) {
            int f = I.fitness;
            if (f == -1) {
                // Compute fitness if not cached
                f = fitness(I.A, d, t, vs, l, non_locating_pairs, max_possible_fitness, is_detecting);
                I.fitness = f; 
            }
            if (f >= max_possible_fitness) {
                return I.A; // Solution found!
            }
            fitnesses.push_back(std::make_pair(f,I));
        }

        // --- Selection ---
        // Sort by fitness
        std::sort(fitnesses.begin(), fitnesses.end(), [](const auto& first, const auto& second) {
            return first.first < second.first;
        });

        // Elitism: Keep the top 50%
        std::vector<Ind_NonRecompute_Fitness> new_vec;
        for (int i=pop_size/2; i < pop_size; i++) {
            new_vec.push_back(fitnesses[i].second);
        }
        pop = new_vec;
        new_vec.clear();

        // --- Crossover & Mutation ---
        // Fill the other 50%
        while (new_vec.size() < pop_size / 2) {
            auto idx1 = any_int(rng) % pop.size();
            auto idx2 = any_int(rng) % pop.size();
            const auto& p1 = pop[idx1];
            const auto& p2 = pop[idx2];

            auto cross_percent = any_int(rng) % 10;
            auto mut_percent = any_int(rng) % 10;
            
            if (cross_percent == 0 && mut_percent < 3) {
                // Crossover + Mutate
                auto new_ind = cross(p1.A,p2.A,d,t,vs,l, rng);
                new_ind = mutate(new_ind,d,t,vs,l, rng);
                Ind_NonRecompute_Fitness true_new_ind;
                true_new_ind.A = new_ind;
                true_new_ind.fitness = -1; // Mark for re-computation
                new_vec.push_back(true_new_ind);
            } else if (cross_percent == 0) {
                // Crossover only
                auto new_ind = cross(p1.A,p2.A,d,t,vs,l, rng);
                Ind_NonRecompute_Fitness true_new_ind;
                true_new_ind.A = new_ind;
                true_new_ind.fitness = -1;
                new_vec.push_back(true_new_ind);
            }
            // (Note: This GA has no "mutate only" path, and a high
            // chance of doing nothing, which is unusual)
        }
        pop.insert(pop.end(), new_vec.begin(), new_vec.end());
    }

    return s; // Failed to find a solution
}


/**
 * @brief Simulated Annealing (SA) search algorithm.
 *
 * Tries to find an array of size 'N' that fixes *all* pairs in
 * 'only_these_pairs'.
 *
 * @param N The number of rows in the array to generate.
 * @param only_these_pairs The list of pairs this array *must* fix.
 * @param rng The thread-local random number generator.
 * @return The successful array, or an empty array on failure.
 */
ca_type try_N_SA(N_type N, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& only_these_pairs, bool is_detecting, std::mt19937& rng) {

    ca_type empty;
    // Start with a random array
    ca_type A = random_array(N, vs.size(), vs);
    
    // SA parameters
    auto temp = 1.0;
    auto rate = 0.99; // Cooling rate
    auto num_iter = 1000;

    // Target fitness: must fix all pairs
    auto required_fitness = only_these_pairs.size();

    auto f = fitness(A, d, t, vs, l, only_these_pairs, required_fitness, is_detecting);
    
    for (int it=0; it<num_iter; it++) {
        
        if (f >= required_fitness) {
            return A; // Solution found!
        }
        
        // Create a new candidate solution by mutating the current one
        auto A_prime = mutate(A, d, t, vs, l, rng);
        auto f_prime = fitness(A_prime, d, t, vs, l, only_these_pairs, required_fitness, is_detecting);
        
        auto diff = f_prime - f;
        
        if (f_prime >= f) {
            // New solution is better, always accept it
            A = A_prime;
            f = f_prime;
        } else {
            // New solution is worse. Accept it with probability e^(-diff/temp)
            auto prob = std::exp(diff / temp); // Note: diff is negative
            if (prob > unif(rng)) { // Use > for prob > random
                A = A_prime;
                f = f_prime;
            }
        }
        // Cool the temperature
        temp = rate * temp;
    }

    if (f >= required_fitness) {
        return A; // Check one last time
    }

    return empty; // Failed to find a solution
}


/**
 * @brief Finds the *minimal* N required to fix a 'percent' of pairs.
 *
 * This function wraps `try_N_SA` and uses an exponential-then-binary
 * search to find the smallest 'N' that works.
 *
 * @param non_locating_pairs The *full* list of pairs.
 * @param percent The *percentage* of pairs from the full list to fix.
 * @return The minimal array that fixes the target subset of pairs.
 */
ca_type go(const d_type& d, const t_type& t, const vs_type& vs, const lambda_type& l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, const double& percent, bool is_detecting, std::mt19937& rng) {
    
    bool succ_first = true;
    ca_type result;

    // Create the subset of pairs to fix in this step
    const std::vector<std::tuple<d_set_type, d_set_type, int>> only_these_pairs(
        non_locating_pairs.begin(), 
        non_locating_pairs.begin() + non_locating_pairs.size() * percent
    );
    
    // Set lower bound for N based on the max separation still needed
    int N = 1;
    for (const auto& [d_set1, d_set2, num] : only_these_pairs) {
        N = std::max(N, l-num);
    }
    if (is_detecting) N = 1; // For detecting, N=1 is always a valid start

    // --- 1. Exponential Search ---
    // Double N until we find *any* solution
    while (true) {
        result = try_N_SA(N, d, t, vs, l, only_these_pairs, is_detecting, rng);
        if (succ_first &&  result.size() > 0) {
            return result; // Succeeded on the first try (N=lower_bound)
        }
        if (result.size() > 0) {
            break; // Found an upper bound
        }
        N *= 2;
        succ_first = false;
    }

    // --- 2. Binary Search ---
    // Now we know a solution exists at 'N', but not at 'N/2'.
    // Binary search between [N/2, N] to find the minimum.
    int N_hi = N;
    int N_lo = N / 2;
    while (N_lo < N_hi) {
        int N_mid = (N_lo + N_hi) / 2;
        auto result2 = try_N_SA(N_mid, d, t, vs, l, only_these_pairs, is_detecting, rng);
        if (result2.size() > 0) {
            // Solution found at N_mid, so this is our new upper bound
            N_hi = N_mid;
            result = result2;
        } else {
            // Failed at N_mid, so the solution must be > N_mid
            N_lo = N_mid + 1;
        }
    }
    return result; // This 'result' is the one for the minimal N (N_hi)
}

/**
 * @brief Checks if 'ind' dominates 'other' (Pareto dominance).
 * An individual dominates if it is no worse in all objectives (N, time)
 * and strictly better in at least one.
 */
bool dominates(const PercentGAFitnessInd& ind, const PercentGAFitnessInd& other) {
    // Note: This implementation is slightly wrong.
    // It should be (ind.N <= other.N && ind.time < other.time) || (ind.N < other.N && ind.time <= other.time)
    // The current version finds the "weakly" non-dominated set.
    return (ind.N <= other.N && 
        ind.time <= other.time);
}

/**
 * @brief Finds the Pareto front from a set of points.
 *
 * @param points The population of solutions (individuals).
 * @return A pair containing:
 * 1. A vector of non-dominated points (the Pareto front).
 * 2. A vector of dominated points.
 */
auto pareto_and_rest(std::vector<PercentGAFitnessInd> points) {
    int candidate_ind_number = 0;
    std::vector<PercentGAFitnessInd> dominated_pts;
    std::vector<PercentGAFitnessInd> pareto;
    
    // This is a simple (but slow, O(n^2)) dominance check
    while (true) {
        auto candidate_ind = points[candidate_ind_number];
        points.erase(points.begin() + candidate_ind_number);
        bool non_dominated = true; 
        int ind_number = 0;
        
        while (points.size() != 0 && ind_number < points.size()) {
            auto ind = points[ind_number];
            if (dominates(candidate_ind, ind)) {
                // Candidate dominates 'ind', so 'ind' is removed
                points.erase(points.begin() + ind_number);
                dominated_pts.push_back(ind);
            } else if (dominates(ind, candidate_ind)) {
                // 'ind' dominates candidate, so candidate is dominated
                non_dominated = false;
                dominated_pts.push_back(candidate_ind);
                ind_number++;
            } else {
                // Neither dominates, move on
                ind_number++;
            }
        }
        if (non_dominated) {
            // Candidate was not dominated by any 'ind'
            pareto.push_back(candidate_ind);
        }
        if (points.size() == 0) {
            break;
        }
    }
    return std::make_pair(pareto, dominated_pts);
}

/**
 * @brief Generates a random "strategy" (individual) for the percent-GA.
 * A strategy is a sorted vector of random percentages, ending in 1.0.
 */
auto generate_rand_percent_individual() {
    std::vector<double> percents;
    int rand_length = ind_size(rng); // 10-30 stages
    for (int i=0; i<rand_length; i++) {
        percents.push_back(unif(rng));
    }
    std::sort(percents.begin(), percents.end());
    percents[0] = 0.001; // Ensure at least 0.1%
    percents.push_back(1.0); // Ensure it always finishes
    PercentGAFitnessInd result;
    result.percents = percents;
    return result;
}

/**
 * @brief Runs the "meta" Genetic Algorithm to find the best strategy.
 *
 * This function is templated to allow either serial (std::execution::seq)
 * or parallel (std::execution::par) execution.
 *
 * @param policy The execution policy (serial or parallel).
 * @param ... Other GA parameters.
 * @return The Pareto front of the *best strategies* found.
 */
template<typename Policy>
auto run_ga_with_policy(
    Policy policy,
    d_type d, 
    t_type t, 
    const vs_type& vs, 
    lambda_type l, 
    const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, 
    bool use_default_percents, 
    bool is_detecting) {

    // --- Option 1: Use a single, hardcoded default strategy ---
    if (use_default_percents) {
        std::mt19937 main_rng(std::random_device{}());
        const std::vector<double> percents = {0.021576,0.021576,0.022644,0.030792,0.090424,0.071014,0.083679,0.172455,0.220123,0.415283,1.000000};
        
        auto non_locating_pairs_copy = non_locating_pairs;
        int num_rows = 0;
        ca_type all_ga_rows; // <<< --- ADD THIS LINE
        auto start = high_resolution_clock::now();
        
        // Execute the strategy step-by-step
        for (const auto& percent : percents) { 
            std::vector<std::tuple<d_set_type, d_set_type, int>> new_non_locating_pairs;

            // Find the minimal array 'ga_rows' to fix this percentage of pairs
            auto ga_rows = go(d,t,vs,l,non_locating_pairs_copy,percent, is_detecting, main_rng);
            num_rows += ga_rows.size();
            all_ga_rows.insert(all_ga_rows.end(), ga_rows.begin(), ga_rows.end());
            
            // --- Update the remaining pairs ---
            // This block checks which pairs were fixed by 'ga_rows' and
            // updates the 'num_times_sep_already' count for the rest.
            std::unordered_map<d_set_type, std::vector<N_type>, DSetHasher> rows_map_for_update;
            auto rows_of_dset_in_ga = [&](const d_set_type& d_set) {
                if (rows_map_for_update.count(d_set)) {
                    return rows_map_for_update.at(d_set);
                }
                robin_hood::unordered_set<N_type> the_rows;
                for (const auto& interaction : d_set) {
                    const auto& rows = rows_of_interaction(interaction, ga_rows);
                    the_rows.insert(rows.begin(), rows.end());
                }
                std::vector<int> vrows(the_rows.begin(), the_rows.end());
                std::sort(vrows.begin(), vrows.end());
                rows_map_for_update[d_set] = vrows;
                return vrows;
            };

            for (const auto& [dset_1, dset_2, num_times_sep_already] : non_locating_pairs_copy) {
                auto rows1 = rows_of_dset_in_ga(dset_1);
                auto rows2 = rows_of_dset_in_ga(dset_2);

                int n = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());
                auto required_separation = is_detecting ? 1 : l;

                // If still not fixed, add to the next generation's list
                if (num_times_sep_already + n < required_separation) {
                    new_non_locating_pairs.push_back({dset_1, dset_2, num_times_sep_already + n});
                }
            }
            non_locating_pairs_copy = new_non_locating_pairs;
            // --- End of update ---

            if (non_locating_pairs_copy.size() == 0) {
                break; // All pairs fixed
            }

            std::cout << "Added " << num_rows << " rows, there are " << non_locating_pairs_copy.size() << " remaining pairs\n";
            new_non_locating_pairs.clear();
        }
        auto stop = high_resolution_clock::now();
        auto total_time = duration_cast<milliseconds>(stop-start).count();
        
        // Return the single result
        std::vector<PercentGAFitnessInd> result;
        PercentGAFitnessInd ind;
        ind.N = num_rows;
        ind.percents = percents;
        ind.time = total_time;
        ind.generated_rows = all_ga_rows;
        result.push_back(ind);
        return result;
    }

    // --- Option 2: Run the full Genetic Algorithm to *find* a good strategy ---
    int pop_size = 100;
    int num_gens = 50;

    std::vector<PercentGAFitnessInd> result;

    // Initialize population with random strategies
    std::vector<PercentGAFitnessInd> pop;
    for (int i=0; i<pop_size; i++) {
        pop.push_back(generate_rand_percent_individual());
    }

    // --- GA Generations Loop ---
    for (int gen=0; gen<num_gens; gen++) {
        std::cout << "Generation #" << gen << "\n";

        // --- Fitness Evaluation (Parallel) ---
        // This loop calculates the fitness (N, time) for every individual (strategy)
        // in the population, using the parallel/sequential policy.
        std::for_each(policy, pop.begin(), pop.end(), 
            [&](PercentGAFitnessInd& I) { // Note: 'I' is one strategy
            
            if (I.N != -1 && I.time != -1) {
                return; // Fitness already known, skip
            }

            // CRITICAL: Each thread needs its *own* private RNG.
            std::mt19937 thread_rng(std::random_device{}());
            
            long long num_rows = 0;
            ca_type all_ga_rows;
            auto non_locating_pairs_copy = non_locating_pairs;
            auto start = high_resolution_clock::now();
            
            // Execute the strategy (I.percents) step-by-step
            for (const auto& percent : I.percents) { 
                std::vector<std::tuple<d_set_type, d_set_type, int>> new_non_locating_pairs;

                // Pass the thread-local RNG to the 'go' function
                auto ga_rows = go(d, t, vs, l, non_locating_pairs_copy, percent, is_detecting, thread_rng);
                num_rows += ga_rows.size();
                all_ga_rows.insert(all_ga_rows.end(), ga_rows.begin(), ga_rows.end());

                // --- Update remaining pairs (logic is identical to the default block) ---
                std::unordered_map<d_set_type, std::vector<N_type>, DSetHasher> rows_map_for_update;
                auto rows_of_dset_in_ga = [&](const d_set_type& d_set) {
                    if (rows_map_for_update.count(d_set)) {
                        return rows_map_for_update.at(d_set);
                    }
                    robin_hood::unordered_set<N_type> the_rows;
                    for (const auto& interaction : d_set) {
                        const auto& rows = rows_of_interaction(interaction, ga_rows);
                        the_rows.insert(rows.begin(), rows.end());
                    }
                    std::vector<int> vrows(the_rows.begin(), the_rows.end());
                    std::sort(vrows.begin(), vrows.end());
                    rows_map_for_update[d_set] = vrows;
                    return vrows;
                };

                for (const auto& [dset_1, dset_2, num_times_sep_already] : non_locating_pairs_copy) {
                    auto rows1 = rows_of_dset_in_ga(dset_1);
                    auto rows2 = rows_of_dset_in_ga(dset_2);
                    int n = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());
                    auto required_separation = is_detecting ? 1 : l;
                    if (num_times_sep_already + n < required_separation) {
                        new_non_locating_pairs.push_back({dset_1, dset_2, num_times_sep_already + n});
                    }
                }
                non_locating_pairs_copy = new_non_locating_pairs;
                // --- End of update ---
            }
            auto stop = high_resolution_clock::now();

            // Update the individual's fitness values.
            I.N = num_rows;
            I.time = duration_cast<milliseconds>(stop-start).count();
            I.generated_rows = all_ga_rows;
        });
        // --- End of Parallel Fitness Evaluation ---

        
        std::vector<PercentGAFitnessInd> fitnesses = pop;

        // --- Selection (Multi-objective) ---
        // Find the Pareto front of the current population
        auto target_size = pop_size/2;
        auto [pareto, rest] = pareto_and_rest(fitnesses);
        
        // Print the current best results
        for (auto& ind : pareto) { // Use 'ind' instead of structured binding
            std::cout << ind.N << "," << ind.time << ",";
            print_vec(ind.percents);
            std::cout << "\n";
        }
        
        // The new population starts with the Pareto front
        std::vector<PercentGAFitnessInd> new_pop(pareto.begin(), pareto.end());
        result = pareto; // Save the best front found so far

        // Fill the rest of the new population with "second-best" fronts
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
            pop.push_back(elem);
        }
        new_pop.clear();

        // --- Crossover & Mutation (on the 'percents' vectors) ---
        std::vector<PercentGAFitnessInd> children;

        // Crossover
        while (children.size() < pop_size/2) {
            auto rand_p1 = pop[any_int(rng) % pop.size()]; // Parent 1 strategy
            auto rand_p2 = pop[any_int(rng) % pop.size()]; // Parent 2 strategy
            
            // One-point crossover on the vector of percentages
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

        // Mutate
        for (auto& child : children) {
            auto r = any_int(rng) % 10;
            if (r == 0) {
                // "split": insert a new random percent
                auto rand_idx = any_int(rng) % (child.percents.size()-1); 
                auto rand_val = unif(rng);
                child.percents.insert(child.percents.begin() + rand_idx, rand_val);
            } else if (r == 1) {
                // "join": remove a percent
                auto rand_idx = any_int(rng) % (child.percents.size()-1); 
                child.percents.erase(child.percents.begin() + rand_idx);
            } else if (r == 2) {
                // "mutate": change one percent
                auto rand_idx = any_int(rng) % (child.percents.size()-1); 
                child.percents[rand_idx] = unif(rng);
            }
            // Add child to population
            pop.push_back(child);
        }
    }
    return result; // Return the last-computed Pareto front
}

/**
 * @brief Public wrapper function for the 'percent_GA'.
 *
 * This function selects the execution policy (parallel or sequential)
 * based on the input string and calls the templated 'run_ga_with_policy'.
 */
std::vector<PercentGAFitnessInd> percent_GA(d_type d, t_type t, const vs_type& vs, const lambda_type& l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, bool use_default_percents, bool is_detecting, const std::string& execution_policy) {
    if (execution_policy == "parallel") {
        std::cout << "Running GA with std::execution::par\n";
        return run_ga_with_policy(std::execution::par, d, t, vs, l, non_locating_pairs, use_default_percents, is_detecting);
    } else {
        std::cout << "Running GA with std::execution::seq\n";
        return run_ga_with_policy(std::execution::seq, d, t, vs, l, non_locating_pairs, use_default_percents, is_detecting);
    }
}