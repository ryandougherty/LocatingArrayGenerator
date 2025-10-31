#pragma once
#include <vector>
#include <cmath>
#include <cstdint>
#include "../combinations.hpp"
#include "../enumerate.hpp"
#include <chrono>
#include <random>
#include <iostream>
#include <fstream>

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

extern std::mt19937_64 rng;//(0); // set seed
extern std::uniform_real_distribution<double> unif;//(0, 1);
extern std::uniform_int_distribution<int> ind_size;//(10, 30);
extern std::uniform_int_distribution<int> any_int;

auto interaction_to_str(const interaction_type& I);
void print_interaction(const interaction_type& I);
void print_d_set(const d_set_type& D);
auto d_set_to_str(const d_set_type& D);
void print_array(const ca_type& A);
void print_array(std::vector<std::vector<int>> &array);
void write_to_file(std::vector<std::vector<int>> &array, std::string file_name);

template <typename T>
void print_vec(const std::vector<T>& vec) {
    for (const auto& elem : vec) {
        std::cout << elem << ", ";
    }
    std::cout << std::endl; // Add a newline for cleaner output
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& p) {
    os << "(";
    // You might need to handle how to print the contents of the pair's vectors
    // For now, let's just indicate their size as an example
    os << "cols: " << p.first.size() << ", vals: " << p.second.size(); 
    os << ")";
    return os;
}

long long comb(unsigned n, unsigned k);
double calc_p(const int t, const int v);
ca_type random_array(const N_type N, const k_type k, const vs_type& vs);

// Templated VectorHasher to work with any vector type
struct VectorHasher {
    template <typename T>
    int operator()(const std::vector<T>& V) const {
        int hash = V.size();
        for(auto const &i : V) {
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

/**
 * @brief Saves a covering array to a CSV file.
 *
 * @param A The array (ca_type) to save.
 * @param filename The name of the file to create (e.g., "output.csv").
 */
void save_array_to_csv(const ca_type& A, const std::string& filename);

/**
 * @brief Calculates the binomial coefficient C(n, k) or "n choose k".
 *
 * Used by the Conditional Expectation algorithm for probability calculations.
 */
namespace math {
	long long combinations(int n, int k);
}