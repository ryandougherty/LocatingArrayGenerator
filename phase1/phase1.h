#pragma once// Header Guard
#include "../utils/utils.h"
#include "../utils/interaction_codec.h"
#include "../robin_hood.h"
#include "../zip.hpp"
#include "../combinations.hpp"
#include "../enumerate.hpp"
#include "../range.hpp"
#include <map>

auto get_interactions(const t_type t, const vs_type& vs, bool t_bar);
bool is_subset(const std::vector<int>& a, const std::vector<int>& b);

// A helper function for explicit and safe comparison of interaction_type objects

std::vector<undist_pair_type> find_non_detecting_sets(const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar, int X, InteractionCodec& codec);

std::vector<undist_pair_type> find_non_locating_sets(const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar, int X, InteractionCodec& codec);

robin_hood::unordered_flat_set<int> rows_of_interaction(const interaction_type& I, const ca_type& A);


template <class InputIterator1, class InputIterator2>
int size_of_symmetric_difference(InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2)
{
    int size = 0;
    while (true)
    {
        if (first1 == last1) return std::distance(first2, last2) + size;
        if (first2 == last2) return std::distance(first1, last1) + size;

        if (*first1 < *first2) { ++first1; ++size; }
        else if (*first2 < *first1) { ++first2; ++size; }
        else { ++first1; ++first2; }
    }
}

/** |A \ B| for sorted ranges — elements in A but not in B. Used for detecting arrays. */
template <class InputIterator1, class InputIterator2>
int size_of_set_minus(InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2)
{
    int size = 0;
    while (first1 != last1)
    {
        if (first2 == last2) return std::distance(first1, last1) + size;
        if (*first1 < *first2) { ++first1; ++size; }
        else if (*first2 < *first1) { ++first2; }
        else { ++first1; ++first2; }
    }
    return size;
}