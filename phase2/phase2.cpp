#include "../utils/utils.h"
#include "phase2.h"
#include "../phase1/phase1.h"
#include <execution>

int fitness(const ca_type& ind, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, const int& threshold, bool is_detecting) {
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
        auto requirement = (is_detecting ? 1 : l) - num_times_sep_already;
        auto rows1 = rows_of_dset(dset_1);
        auto rows2 = rows_of_dset(dset_2);
        int n = size_of_symmetric_difference(rows1.begin(), rows1.end(), rows2.begin(), rows2.end());
        
        if (n >= requirement) {
            score += 1;
        }

        if (score >= threshold) {
            return threshold + 1;
        }
    }
    return score;
}

ca_type cross(const ca_type& p1, const ca_type& p2, d_type d, t_type t, const vs_type& vs, lambda_type l, std::mt19937& rng) {
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

ca_type mutate(const ca_type& p1, d_type d, t_type t, const vs_type& vs, lambda_type l, std::mt19937& rng) {
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

ca_type try_N(N_type N, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, double percent, bool is_detecting, std::mt19937& rng) {

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
                f = fitness(I.A, d, t, vs, l, non_locating_pairs, max_possible_fitness, is_detecting);
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
                auto new_ind = cross(p1.A,p2.A,d,t,vs,l, rng);
                new_ind = mutate(new_ind,d,t,vs,l, rng);
                Ind_NonRecompute_Fitness true_new_ind;
                true_new_ind.A = new_ind;
                true_new_ind.fitness = -1;
                new_vec.push_back(true_new_ind);
            } else if (cross_percent == 0) {
                auto new_ind = cross(p1.A,p2.A,d,t,vs,l, rng);
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




ca_type try_N_SA(N_type N, d_type d, t_type t, const vs_type& vs, lambda_type l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& only_these_pairs, bool is_detecting, std::mt19937& rng) {



    ca_type empty;
    ca_type A = random_array(N, vs.size(), vs);
    auto temp = 1.0;
    auto rate = 0.99;
    auto num_iter = 1000;

    auto required_fitness = only_these_pairs.size();

    auto f = fitness(A, d, t, vs, l, only_these_pairs, required_fitness, is_detecting);
    for (int it=0; it<num_iter; it++) {
        
        if (f >= required_fitness) {
            return A;
        }
        auto A_prime = mutate(A, d, t, vs, l, rng);
        auto f_prime = fitness(A_prime, d, t, vs, l, only_these_pairs, required_fitness, is_detecting);
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

ca_type go(const d_type& d, const t_type& t, const vs_type& vs, const lambda_type& l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, const double& percent, bool is_detecting, std::mt19937& rng) {
    
    
    bool succ_first = true;
    ca_type result;
    std::size_t howManyPairs;

    const std::vector<std::tuple<d_set_type, d_set_type, int>> only_these_pairs(non_locating_pairs.begin(), non_locating_pairs.begin() + non_locating_pairs.size() * percent);
    
    int N = 1;
    for (const auto& [d_set1, d_set2, num] : only_these_pairs) {
        N = std::max(N, l-num);
    }

    while (true) {

        result = try_N_SA(N, d, t, vs, l, only_these_pairs, is_detecting, rng);
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
        auto result2 = try_N_SA(N_mid, d, t, vs, l, only_these_pairs, is_detecting, rng);
        if (result2.size() > 0) {
            N_hi = N_mid;
            result = result2;
        } else {
            N_lo = N_mid + 1;
        }
    }
    return result;
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

// std::vector<PercentGAFitnessInd> percent_GA(d_type d, t_type t, const vs_type& vs, const lambda_type& l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, bool use_default_percents, bool is_detecting, std:string policy) {
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

    if (use_default_percents) {
        std::mt19937 main_rng(std::random_device{}());
        const std::vector<double> percents = {0.021576,0.021576,0.022644,0.030792,0.090424,0.071014,0.083679,0.172455,0.220123,0.415283,1.000000};
//{0.001000,0.002625,0.010040,0.017456,0.631592,0.094666,0.151855,0.167999,0.172241,1.000000};
        auto non_locating_pairs_copy = non_locating_pairs;
        int num_rows = 0;
        auto start = high_resolution_clock::now();
        for (const auto& percent : percents) { 
        // std::for_each(std::execution::par, percents.begin(), percents.end(), 
            // [&](const auto& percent) {
            std::vector<std::tuple<d_set_type, d_set_type, int>> new_non_locating_pairs;

            auto ga_rows = go(d,t,vs,l,non_locating_pairs_copy,percent, is_detecting, main_rng);
            num_rows += ga_rows.size();
            // *** START OF THE FIX ***
            // This block now correctly updates the list of remaining pairs for both modes.
            
            // Define a helper lambda to get rows for a d-set within the generated array
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

            if (non_locating_pairs_copy.size() == 0) {
                break;
            }

            std::cout << "Added " << num_rows << " rows, there are " << non_locating_pairs_copy.size() << " remaining pairs\n";
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

        // This is the main change: We use std::for_each with a parallel policy
        // to calculate the fitness for every individual in the population concurrently.
        std::for_each(policy, pop.begin(), pop.end(), 
            [&](PercentGAFitnessInd& I) {
            
            // If fitness is already known, we can skip this individual.
            if (I.N != -1 && I.time != -1) {
                return;
            }

            // CRITICAL: Each thread must have its own private random number generator
            // to prevent data races. We seed it with a true hardware random device.
            std::mt19937 thread_rng(std::random_device{}());
            
            long long num_rows = 0;
            auto non_locating_pairs_copy = non_locating_pairs;

            auto start = high_resolution_clock::now();
            for (const auto& percent : I.percents) { 
                std::vector<std::tuple<d_set_type, d_set_type, int>> new_non_locating_pairs;

                // Pass the thread-local RNG down the call stack.
                auto ga_rows = go(d, t, vs, l, non_locating_pairs_copy, percent, is_detecting, thread_rng);
                num_rows += ga_rows.size();

                // This logic is safe because all variables are local to this thread's execution.
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
            }
            auto stop = high_resolution_clock::now();

            // Update the individual's computed fitness values. This is thread-safe
            // because each thread is writing to a different individual 'I'.
            I.N = num_rows;
            I.time = duration_cast<milliseconds>(stop-start).count();
        });

        // After the parallel loop finishes, 'pop' is fully updated.
        // We now copy it to 'fitnesses' to proceed with selection.
        std::vector<PercentGAFitnessInd> fitnesses = pop;

        // The rest of your genetic algorithm (selection, crossover, mutation) remains the same.
        auto target_size = pop_size/2;
        auto [pareto, rest] = pareto_and_rest(fitnesses);
        for (auto& [percents, N, time] : pareto) {
            std::cout << N << "," << time << ",";
            print_vec(percents);
            std::cout << "\n";
        }
        std::vector<PercentGAFitnessInd> new_pop(pareto.begin(), pareto.end());

        result = pareto;

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

std::vector<PercentGAFitnessInd> percent_GA(d_type d, t_type t, const vs_type& vs, const lambda_type& l, const std::vector<std::tuple<d_set_type, d_set_type, int>>& non_locating_pairs, bool use_default_percents, bool is_detecting, const std::string& execution_policy) {
    if (execution_policy == "parallel") {
        return run_ga_with_policy(std::execution::par, d, t, vs, l, non_locating_pairs, use_default_percents, is_detecting);
    } else {
        return run_ga_with_policy(std::execution::seq, d, t, vs, l, non_locating_pairs, use_default_percents, is_detecting);
    }
}