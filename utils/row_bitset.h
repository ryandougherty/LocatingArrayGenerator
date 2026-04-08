/**
 * row_bitset.h
 *
 * Lightweight dynamic bitset for representing row sets.
 * Bit i is set if row i covers the interaction/d-set.
 *
 * Key operations (N = number of rows, W = ceil(N/64) words):
 *   - symmetric_diff_popcount: W XORs + W popcounts  (vs O(N) merge scan)
 *   - union_with:              W ORs                  (vs O(N) set insert)
 *   - set/clear/test:          O(1)
 *
 * For N=200, W=4 — a symmetric difference is 4 XORs + 4 popcounts.
 */

#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>

struct RowBitset {
    std::vector<uint64_t> words;
    int num_rows = 0;

    RowBitset() = default;
    explicit RowBitset(int n) : words(((unsigned)n + 63u) / 64u, 0ULL), num_rows(n) {}

    void set(int i)          { words[(unsigned)i >> 6] |=  (1ULL << (i & 63)); }
    void clear(int i)        { words[(unsigned)i >> 6] &= ~(1ULL << (i & 63)); }
    bool test(int i) const   { return (words[(unsigned)i >> 6] >> (i & 63)) & 1; }
    void reset()             { std::fill(words.begin(), words.end(), 0ULL); }

    int popcount() const {
        int c = 0;
        for (auto w : words) c += __builtin_popcountll(w);
        return c;
    }

    /** |A XOR B| — number of rows in exactly one of the two sets (locating). */
    int symmetric_diff_popcount(const RowBitset& o) const {
        int c = 0;
        for (size_t i = 0; i < words.size(); ++i)
            c += __builtin_popcountll(words[i] ^ o.words[i]);
        return c;
    }

    /** |A \ B| — number of rows in A but not B (detecting). */
    int andnot_popcount(const RowBitset& o) const {
        int c = 0;
        for (size_t i = 0; i < words.size(); ++i)
            c += __builtin_popcountll(words[i] & ~o.words[i]);
        return c;
    }

    /**
     * Unified separation metric.
     *   Locating:  |A △ B| = popcount(A XOR B)
     *   Detecting: |A \ B| = popcount(A AND NOT B)
     *     where A = rows of interaction X, B = rows of d-set T
     */
    int separation(const RowBitset& o, bool detecting) const {
        if (detecting) return andnot_popcount(o);
        else           return symmetric_diff_popcount(o);
    }

    /** this |= other */
    void union_with(const RowBitset& o) {
        for (size_t i = 0; i < words.size(); ++i)
            words[i] |= o.words[i];
    }

    /** Build from OR of several other bitsets. */
    void build_union(const RowBitset* srcs, size_t count) {
        reset();
        for (size_t s = 0; s < count; ++s)
            union_with(srcs[s]);
    }

    bool operator==(const RowBitset& o) const { return words == o.words; }
    bool operator!=(const RowBitset& o) const { return words != o.words; }
};
