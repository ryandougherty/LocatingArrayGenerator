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
#include "utils/utils.h"
#include "phase1/phase1.h"
#include "phase2/phase2.h"

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
    //other papers with these for comparison***
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

    if (argc != 4) {
        std::cerr << "Usage: ./LocAG <name of config> <array_type> <execution_policy>\n";
        std::cerr << "array_type can be 'locating' or 'detecting'\n";
        return -1;
    }

    std::string array_type = argv[2];
    const std::string policy = argv[3];
    

    for (d_type d = 1; d <= 1; d++) {
        for (t_type t = 2; t <= 2; t++) {
            for (lambda_type lambda = 1; lambda <= 4; lambda++) {
                const std::string config_name = argv[1];
                const auto& [vs, filename] = lookup_config_and_params(config_name, t, lambda);

                const bool d_bar = true;
                const bool t_bar = true;
                bool is_detecting = false;

                assert(d < *std::min_element(vs.begin(), vs.end()));

                std::cout << "------------d=" << std::to_string(d) << ", t=" << std::to_string(t) << ", lambda=" << std::to_string(lambda) << ", filename=" << filename << "------------\n";

                auto start = high_resolution_clock::now();
                ca_type A = read_ca_from_cagen(filename, vs);
                // std::cout << "Read file with " << A.size() << " rows.\n";

                // Finds initial non_locating_pairs
                std::vector<std::tuple<d_set_type, d_set_type, int>> non_valid_pairs;
                if (array_type == "locating") {
                    non_valid_pairs = find_non_locating_sets(A, t, vs, lambda, d, d_bar, t_bar);
                }
                else if (array_type == "detecting") {
                    is_detecting = true;
                    non_valid_pairs = find_non_detecting_sets(A, t, vs, lambda, d, d_bar, t_bar);
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

                /* ------------------------------------- Stage 2: GA ------------------------------ */

                if (non_valid_pairs.size() == 0) {
                    break;
                }

                auto pareto = percent_GA(d,t,vs,lambda,non_valid_pairs,true, is_detecting, policy);
                for (const auto& [percents, num_rows, time] : pareto) {
                    std::cout << "N total=" << first_stage_N + num_rows << ", Time total=" << first_stage_time + time << ", percents=";
                    print_vec(percents);
                    std::cout << "\n";
                }
                
            }
        }
    }
    

}

//Do checks in parallel and for main loop
//Do GA double for loop in parallel
//Experimenting (look at papers on CA and software testing) look in transactions in sofwater engineering for research questions
//Research Questions: 
//1. Is this tool faster than previous tools? Is it more memory efficient? Does it run better in parallel? 
//Start writing a summary of the method (based on other paper) and new work