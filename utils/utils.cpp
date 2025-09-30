#include "utils.h"

std::mt19937_64 rng(0); // set seed
std::uniform_real_distribution<double> unif(0, 1);
std::uniform_int_distribution<int> ind_size(10, 30);
std::uniform_int_distribution<int> any_int;

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

void print_interaction(const interaction_type& I) {
    std::cout << interaction_to_str(I);
}

// template <typename T>
// void print_vec(const std::vector<T>& vec) {
//     for (const auto& elem : vec) {
//         std::cout << elem << ", ";
//     }
//     std::cout << std::endl; // Add a newline for cleaner output
// }

void print_d_set(const d_set_type& D) {
    for (const auto& I : D) {
        print_interaction(I);
    }
}

auto d_set_to_str(const d_set_type& D) {
    std::string result;
    for (const auto& I : D) {
        result += interaction_to_str(I) + " ";
    }
    return result;
}

void print_array(const ca_type& A) {
    for (const auto& row : A) {
        for (const auto& value : row) {
            std::cout << std::to_string(value) << ' ';
        }
        std::cout << '\n';
    }
}


long long comb(unsigned n, unsigned k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    long long result = n;
    for (int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

double calc_p(const int t, const int v) {
    return 1 / pow(v, t);
}

ca_type random_array(const N_type N, const k_type k, const vs_type& vs) {
    ca_type to_return(N, std::vector<v_type>(k, 0));
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < k; col++) {
            std::uniform_int_distribution<int> distribution(0, vs[col] - 1);
            to_return[row][col] = distribution(rng);
        }
    }
    return to_return;
}