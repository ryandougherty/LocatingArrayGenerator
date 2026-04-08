/**
 * interaction_codec.h
 *
 * Two-level combinatorial encoding for (at most d)-sets of interactions,
 * adapted for locating and detecting arrays with mixed levels (vs_type).
 *
 * Level 1: interaction_type <-> interaction_id (uint64_t)
 *   Built from the enumeration produced by get_interactions(). Each
 *   interaction gets a sequential index 0, 1, ..., T-1.
 *
 * Level 2: d_set (sorted vector of interaction_ids) <-> d_set_id (uint64_t)
 *   Uses the layered combinatorial number system (colexicographic).
 *   rank(S) = layer_offset[|S|] + colex_rank(S)
 *   where layer_offset[j] = sum_{i=0}^{j-1} C(T, i).
 */

#pragma once

#include "utils.h" // for interaction_type, d_set_type, k_type, v_type, vs_type, t_type, d_type
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <iostream>

using interaction_id = uint64_t;
using d_set_id = uint64_t;

// A pair of d-set IDs + separation count — replaces tuple<d_set_type, d_set_type, int>
using undist_pair_type = std::tuple<d_set_id, d_set_id, int>;

/**
 * @brief Computes binomial coefficient C(n, k) for uint64_t values.
 *
 * Exact for small k (which is all we need — k <= d, typically <= 3).
 * Performs iterative multiplication/division to avoid overflow where possible.
 */
inline uint64_t binom_u64(uint64_t n, int k) {
    if (k < 0 || (uint64_t)k > n) return 0;
    if (k == 0) return 1;
    if ((uint64_t)k > n - (uint64_t)k) k = (int)(n - (uint64_t)k); // symmetry
    uint64_t result = 1;
    for (int i = 0; i < k; ++i) {
        result = result * (n - (uint64_t)i) / (uint64_t)(i + 1);
    }
    return result;
}


/**
 * @brief Encodes/decodes interactions and d-sets as integers.
 *
 * After calling init(), the codec provides:
 *   - encode_interaction(interaction_type) -> interaction_id
 *   - decode_interaction(interaction_id) -> interaction_type
 *   - encode_d_set(d_set_type) -> d_set_id
 *   - decode_d_set(d_set_id) -> d_set_type
 *   - encode_d_set_ids(sorted vector<interaction_id>) -> d_set_id
 *   - decode_d_set_ids(d_set_id) -> sorted vector<interaction_id>
 */
struct InteractionCodec {

    // ----- Level 1 data: interaction <-> interaction_id -----
    std::vector<interaction_type> id_to_interaction;
    std::unordered_map<interaction_type, interaction_id, InteractionHasher> interaction_to_id_map;
    uint64_t T = 0; // total number of interactions

    // ----- Level 2 data: d-set encoding parameters -----
    int d_max = 0;                      // max d-set size
    std::vector<uint64_t> layer_offset; // layer_offset[j] = sum_{i=0}^{j-1} C(T, i)
    uint64_t M = 0;                     // total number of encodable d-sets

    /**
     * @brief Initialize the codec from a pre-generated list of interactions.
     *
     * @param interactions All interactions (from get_interactions()).
     *                     The order defines the integer encoding.
     * @param d            Maximum d-set size.
     */
    void init(const std::vector<interaction_type>& interactions, int d) {
        d_max = d;
        T = interactions.size();

        // Build bidirectional maps
        id_to_interaction = interactions;
        interaction_to_id_map.clear();
        interaction_to_id_map.reserve(T);
        for (uint64_t i = 0; i < T; ++i) {
            interaction_to_id_map[interactions[i]] = i;
        }

        // Precompute layer offsets for Level 2
        layer_offset.resize(d_max + 2, 0);
        for (int j = 0; j <= d_max; ++j) {
            layer_offset[j + 1] = layer_offset[j] + binom_u64(T, j);
        }
        M = layer_offset[d_max + 1];
    }

    // ===================== Level 1: Single interaction =====================

    interaction_id encode_interaction(const interaction_type& ix) const {
        auto it = interaction_to_id_map.find(ix);
        assert(it != interaction_to_id_map.end() && "Interaction not found in codec");
        return it->second;
    }

    const interaction_type& decode_interaction(interaction_id id) const {
        assert(id < T && "Interaction ID out of range");
        return id_to_interaction[id];
    }

    // ===================== Level 2: d-set of interactions ==================

    /**
     * @brief Encode a d_set_type (vector of interaction_type) to a d_set_id.
     */
    d_set_id encode_d_set(const d_set_type& dset) const {
        std::vector<interaction_id> ids(dset.size());
        for (size_t i = 0; i < dset.size(); ++i) {
            ids[i] = encode_interaction(dset[i]);
        }
        std::sort(ids.begin(), ids.end());
        return encode_d_set_ids(ids);
    }

    /**
     * @brief Encode a sorted vector of interaction_ids to a d_set_id.
     *
     * Uses layered colexicographic ranking:
     *   rank = layer_offset[j] + sum_{i=0}^{j-1} C(ids[i], i+1)
     */
    d_set_id encode_d_set_ids(const std::vector<interaction_id>& sorted_ids) const {
        int j = (int)sorted_ids.size();
        assert(j <= d_max);
        uint64_t r = layer_offset[j];
        for (int i = 0; i < j; ++i) {
            r += binom_u64(sorted_ids[i], i + 1);
        }
        return r;
    }

    /**
     * @brief Decode a d_set_id back to a sorted vector of interaction_ids.
     */
    std::vector<interaction_id> decode_d_set_ids(d_set_id r) const {
        // 1. Determine layer j (the d-set size)
        int j = 0;
        while (j <= d_max && layer_offset[j + 1] <= r) ++j;
        assert(j <= d_max);

        // 2. Subtract layer offset
        r -= layer_offset[j];

        // 3. Colex unrank as j-subset of {0, ..., T-1}
        std::vector<interaction_id> ids(j);
        for (int i = j; i >= 1; --i) {
            // Binary search: find largest s such that C(s, i) <= r
            uint64_t lo = (uint64_t)(i - 1), hi = T - 1;
            while (lo < hi) {
                uint64_t mid = lo + (hi - lo + 1) / 2;
                if (binom_u64(mid, i) <= r)
                    lo = mid;
                else
                    hi = mid - 1;
            }
            ids[i - 1] = lo;
            r -= binom_u64(lo, i);
        }
        return ids;
    }

    /**
     * @brief Decode a d_set_id back to a d_set_type (vector of interaction_type).
     */
    d_set_type decode_d_set(d_set_id id) const {
        auto ids = decode_d_set_ids(id);
        d_set_type result(ids.size());
        for (size_t i = 0; i < ids.size(); ++i) {
            result[i] = decode_interaction(ids[i]);
        }
        return result;
    }
};
