/* ----------------------------------------------------------------------------
 * utils.cpp
 *
 * Provides various utility functions, helper classes, and global objects
 * (like the random number generator) used across the project.
 * ----------------------------------------------------------------------------
 */

#include "utils.h"
#include <fstream>
#include <iostream> // For std::cerr
#include <limits>

// --- Global Random Number Generation ---
// A single global RNG, seeded with 0 for reproducible results.
std::mt19937_64 rng(0); 
// Global distribution for doubles between [0, 1]
std::uniform_real_distribution<double> unif(0, 1);
// Global distribution for GA individual size (10 to 30)
std::uniform_int_distribution<int> ind_size(10, 30);
// Global distribution for any integer (default range)
std::uniform_int_distribution<int> any_int;


/**
 * @brief Converts a single interaction to a string.
 * @param I The interaction to convert.
 * @return A string representation, e.g., "(0,1,),(1,0,)"
 */
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

/**
 * @brief Prints an interaction to std::cout.
 */
void print_interaction(const interaction_type& I) {
    std::cout << interaction_to_str(I);
}

// (Generic print_vec is commented out in the original file)
// template <typename T>
// void print_vec(const std::vector<T>& vec) {
//     for (const auto& elem : vec) {
//         std::cout << elem << ", ";
//     }
//     std::cout << std::endl; // Add a newline for cleaner output
// }

/**
 * @brief Prints a d-set (a vector of interactions) to std::cout.
 */
void print_d_set(const d_set_type& D) {
    for (const auto& I : D) {
        print_interaction(I);
    }
}

/**
 * @brief Converts a d-set to a string.
 */
auto d_set_to_str(const d_set_type& D) {
    std::string result;
    for (const auto& I : D) {
        result += interaction_to_str(I) + " ";
    }
    return result;
}

/**
 * @brief Prints a 2D covering array (ca_type) to std::cout.
 */
void print_array(const ca_type& A) {
    for (const auto& row : A) {
        for (const auto& value : row) {
            std::cout << std::to_string(value) << ' ';
        }
        std::cout << '\n';
    }
}

/**
 * @brief Calculates "n choose k" (combinations).
 *
 * @param n The total number of items.
 * @param k The number of items to choose.
 * @return The number of combinations (nCr).
 */
long long comb(unsigned n, unsigned k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n - k; // Take advantage of C(n,k) == C(n, n-k)
    if (k == 0) return 1;

    long long result = n;
    for (int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

/**
 * @brief Calculates the probability of a t-way interaction.
 * Assumes a uniform random distribution of values.
 *
 * @param t The strength.
 * @param v The number of levels (assumed uniform).
 * @return 1 / (v^t)
 */
double calc_p(const int t, const int v) {
    return 1 / pow(v, t);
}

/**
 * @brief Generates a random N x k array.
 *
 * Each cell (row, col) is filled with a random integer from
 * [0, vs[col]-1], respecting the number of levels for that column.
 *
 * @param N The number of rows.
 * @param k The number of columns.
 * @param vs The vector of levels for each column.
 * @return A randomly generated ca_type.
 */
ca_type random_array(const N_type N, const k_type k, const vs_type& vs) {
    // Initialize an array of N rows and k columns
    ca_type to_return(N, std::vector<v_type>(k, 0));
    
    // Fill each cell
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < k; col++) {
            // Create a distribution for this specific column's levels
            std::uniform_int_distribution<int> distribution(0, vs[col] - 1);
            // Assign a random value
            to_return[row][col] = distribution(rng);
        }
    }
    return to_return;
}

/**
 * @brief Saves a covering array to a CSV file.
 *
 * Each row is a line, with values separated by commas.
 *
 * @param A The array (ca_type) to save.
 * @param filename The name of the file to create (e.g., "output.csv").
 */
void save_array_to_csv(const ca_type& A, const std::string& filename) {
    // Open the output file
    std::ofstream output_file(filename);
    
    // Check if the file opened successfully
    if (!output_file.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
        return;
    }

    if (A.empty()) {
        output_file.close();
        return; // Nothing to write
    }

    // Iterate over each row
    for (size_t r = 0; r < A.size(); ++r) {
        // Iterate over each column in the row
        for (size_t c = 0; c < A[r].size(); ++c) {
            output_file << std::to_string(A[r][c]);
            // Add a comma unless it's the last element
            if (c < A[r].size() - 1) {
                output_file << ",";
            }
        }
        // Add a newline at the end of the row
        output_file << "\n";
    }

    // Close the file
    output_file.close();
}

void print_array(std::vector<std::vector<int>> &array) {
	for (unsigned int i = 0; i < array.size(); i++) {
		for (unsigned int j = 0; j < array[i].size(); j++) {
			std::cout << array[i][j];
		}
		std::cout << std::endl;
	}
}

void write_to_file(std::vector<std::vector<int>> &array, std::string file_name) {
	std::ofstream output_file(file_name);
	for (unsigned int i = 0; i < array.size(); i++) {
		for (unsigned int j = 0; j < array[i].size(); j++) {
			output_file << array[i][j];
		}
		output_file << std::endl;
	}
	output_file.close();
}

/**
 * @brief Calculates the binomial coefficient C(n, k) or "n choose k".
 */
namespace math {
    long long combinations(int n, int k) {
        if (k < 0 || k > n) {
            return 0;
        }
        if (k == 0 || k == n) {
            return 1;
        }
        // Take advantage of symmetry C(n, k) = C(n, n-k)
        if (k > n / 2) {
            k = n - k;
        }
        
        long long res = 1;
        for (int i = 1; i <= k; ++i) {
            if (res > std::numeric_limits<long long>::max() / (n - i + 1)) {
                // Handle potential overflow
                std::cerr << "Warning: Overflow in combinations(" << n << ", " << k << ")" << std::endl;
                return std::numeric_limits<long long>::max();
            }
            res = res * (n - i + 1) / i;
        }
        return res;
    }
}