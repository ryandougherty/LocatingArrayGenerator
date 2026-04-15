/* ----------------------------------------------------------------------------
 * phase1.cpp
 *
 * Implements Phase 1 of the LocAG algorithm.
 *
 * The primary goal of this phase is to analyze an initial covering array (A)
 * and identify all pairs of "d-sets" (sets of interactions) that
 * violate the locating or detecting property.
 *
 * This is a major bottleneck. To speed this up, it uses a heuristic:
 * 1.  It divides the array 'A' into 'X' horizontal partitions.
 * 2.  For each interaction 'I', it creates a "signature" vector of counts,
 * recording how many times 'I' appears in each partition.
 * 3.  It groups d-sets based on these signatures.
 * 4.  It assumes that two d-sets (D1, D2) can only be non-located/detected
 * if they have *identical* (or very similar) signatures.
 * 5.  This drastically reduces the number of pairs (D1, D2) that require the
 * expensive full symmetric difference check.
 * ----------------------------------------------------------------------------
 */

#include "phase1.h"
#include "../utils/row_bitset.h"

// A vector of counts, one for each partition.
using partitioned_counts = std::vector<int>;

/**
 * @brief Checks if interaction I1 subsumes interaction I2.
 *
 * I1 subsumes I2 iff:
 *   - cols(I1) ⊆ cols(I2), AND
 *   - for every column in I1, the value matches the corresponding value in I2.
 *
 * When I1 subsumes I2, every row covering I2 also covers I1, so R(I2) ⊆ R(I1).
 * This means having both in a d-set is redundant: R({I1, I2}) = R(I1).
 *
 * Example: I1 = (col0=0) subsumes I2 = (col0=0, col1=1),
 * because any row with col0=0 AND col1=1 also has col0=0.
 */
static bool interaction_subsumes(const interaction_type& I1, const interaction_type& I2) {
    // I1 can only subsume I2 if I1 has fewer or equal columns
    if (I1.first.size() > I2.first.size()) return false;

    // For every (col, val) in I1, check that I2 has the same (col, val)
    for (size_t a = 0; a < I1.first.size(); ++a) {
        bool found = false;
        for (size_t b = 0; b < I2.first.size(); ++b) {
            if (I1.first[a] == I2.first[b]) {
                if (I1.second[a] == I2.second[b]) {
                    found = true;
                }
                break; // column matched, value either matched or didn't
            }
        }
        if (!found) return false;
    }
    return true;
}

/**
 * @brief Checks if a d-set contains a redundant interaction.
 *
 * A d-set has a redundant interaction if any interaction I_a subsumes
 * another interaction I_b (where a != b).  In that case, R(I_b) ⊆ R(I_a),
 * so R({I_a, I_b, ...}) = R({I_a, ...}) — the d-set is equivalent to a
 * smaller one that's already enumerated.
 *
 * This is the key filter for t_bar=true with d > 1: it prevents generating
 * d-sets like {(col0=0), (col0=0, col1=1)} whose row sets are identical
 * to the size-1 d-set {(col0=0)}, which creates structurally indistinguishable
 * pairs that Phase 2 can never resolve.
 */
static bool has_redundant_interaction(const d_set_type& d_set) {
    for (size_t a = 0; a < d_set.size(); ++a) {
        for (size_t b = 0; b < d_set.size(); ++b) {
            if (a == b) continue;
            // If I_a subsumes I_b, then I_b is redundant
            if (interaction_subsumes(d_set[a], d_set[b])) {
                return true;
            }
        }
    }
    return false;
}

/**
 * @brief Calculates the "partitioned count signature" for a single interaction.
 *
 * Divides the array 'A' into 'X' horizontal partitions and counts the
 * occurrences of interaction 'I' within each partition.
 *
 * @param I The interaction (cols, vals) to find.
 * @param A The covering array.
 * @param X The number of partitions.
 * @return A vector of size 'X' where each element is the count for that partition.
 */
auto get_partitioned_interaction_counts(const interaction_type& I, const ca_type& A, int X) {
    partitioned_counts counts(X, 0);
    int partition_size = A.size() / X;

    for (const auto& [idx, row] : enumerate(A)) {
        bool match = true;
        // Check if the row matches the interaction
        for (const auto& [col, val] : zip(I.first, I.second)) {
            if (row[col] != val) {
                match = false;
                break;
            }
        }
        if (match) {
            // If it matches, find which partition this row index belongs to
            int partition_index = std::min(X - 1, static_cast<int>(idx / partition_size));
            counts[partition_index]++;
        }
    }
    return counts;
}

/**
 * @brief Finds all row indices in 'A' where interaction 'I' appears.
 *
 * @param I The interaction (cols, vals) to find.
 * @param A The covering array to search.
 * @return A set of row indices.
 */
robin_hood::unordered_flat_set<int> rows_of_interaction(const interaction_type& I, const ca_type& A) {
    const auto& cols = I.first;
    const auto& vals = I.second;
    robin_hood::unordered_flat_set<int> rows_I_appears;
    
    // Iterate over all rows with their index
    for (const auto& [idx, row] : enumerate(A)) {
        bool flag = false; // Becomes true if this row *doesn't* match
        // Check if the row matches the interaction
        for (const auto& [col, val] : zip(cols, vals)) {
            if (row[col] != val) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            // No mismatch found, so this row covers the interaction
            rows_I_appears.insert(idx);
        }
    }
    return rows_I_appears;
}


/**
 * @brief Generates all interactions of strength up to 't'.
 *
 * An interaction is a (column_set, value_tuple) pair.
 * If 't_bar' is true, it generates interactions for strengths 1, 2, ..., t.
 * If 't_bar' is false, it only generates interactions for strength 't'.
 *
 * @param t The maximum strength.
 * @param vs The vector of levels for each column.
 * @param t_bar Whether to include interactions of strength < t.
 * @return A vector of all possible interactions.
 */
auto get_interactions(const t_type t, const vs_type& vs, bool t_bar) {
    // --- Generate all column sets ---
    auto lb = t; // lower bound
    if (t_bar) {
        lb = 1; // If t_bar, start from strength 1
    }
    std::vector<std::vector<k_type>> col_sets;
    for (int i=lb; i<=t; i++) {
        // Get all combinations of columns of size 'i'
        auto cols = combinations(range(vs.size()), i);
        for (const auto& col_set : cols) {
            std::vector<k_type> to_add;
            for (const auto& new_col : col_set) {
                to_add.push_back(new_col);
            }
            col_sets.push_back(to_add);
        }
    }

    // --- Generate all value-tuples for each column set ---
    // This is a (sub-optimal) hardcoded loop for t=1 up to t=6.
    // A more general solution would use iterators or recursion.
    std::vector<interaction_type> interactions;
    for (const auto& col : col_sets) {
        if (col.size() == 1) {
            for (v_type i = 0; i < vs[col[0]]; i++) {
                std::vector<v_type> s{i};
                interaction_type I = std::make_pair(col, s);
                interactions.push_back(I);
            }
        }
        else if (col.size() == 2) {
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    std::vector<v_type> s{i, j};
                    interaction_type I = std::make_pair(col, s);
                    interactions.push_back(I);
                }
            }
        }
        else if (col.size() == 3) {
            // ... (loops for t=3) ...
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    for(v_type k = 0; k < vs[col[2]]; k++) {
                        std::vector<v_type> s{i, j, k};
                        interaction_type I = std::make_pair(col, s);
                        interactions.push_back(I);
                    }
                }
            }
        }
        else if (col.size() == 4) {
            // ... (loops for t=4) ...
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    for(v_type k = 0; k < vs[col[2]]; k++) {
                        for(v_type l = 0; l < vs[col[3]]; l++) {
                            std::vector<v_type> s{i, j, k, l};
                            interaction_type I = std::make_pair(col, s);
                            interactions.push_back(I);
                        }
                    }
                }
            }
        }
        else if (col.size() == 5) {
            // ... (loops for t=5) ...
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    for(v_type k = 0; k < vs[col[2]]; k++) {
                        for(v_type l = 0; l < vs[col[3]]; l++) {
                            for(v_type m = 0; m < vs[col[4]]; m++) {
                                std::vector<v_type> s{i, j, k, l, m};
                                interaction_type I = std::make_pair(col, s);
                                interactions.push_back(I);
                            }
                        }
                    }
                }
            }
        }
        else if (col.size() == 6) {
            // ... (loops for t=6) ...
            for (v_type i = 0; i < vs[col[0]]; i++) {
                for (v_type j = 0; j < vs[col[1]]; j++) {
                    for(v_type k = 0; k < vs[col[2]]; k++) {
                        for(v_type l = 0; l < vs[col[3]]; l++) {
                            for(v_type m = 0; m < vs[col[4]]; m++) {
                                for (v_type n = 0; n < vs[col[5]]; n++) {
                                    std::vector<v_type> s{i, j, k, l, m, n};
                                    interaction_type I = std::make_pair(col, s);
                                    interactions.push_back(I);
                                }
                            }
                        }
                    }
                }
            }
        } 
        else {
            std::cerr << "Error, Invalid t (t > 6 is not supported by get_interactions)\n";
            abort();
        }
    }
    return interactions;
}

/**
 * @brief Checks if vector 'a' is a subset of vector 'b'.
 * Assumes both vectors are sorted.
 */
bool is_subset(const std::vector<int>& a, const std::vector<int>& b) {
    return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

/**
 * @brief Lexicographical comparison for two interactions.
 * Used for sorting or as a key in std::map.
 */
bool compare_interactions(const interaction_type& a, const interaction_type& b) {
    // First, compare the size of the column vectors
    if (a.first.size() < b.first.size()) return true;
    if (a.first.size() > b.first.size()) return false;

    // If sizes are the same, compare the column vector elements
    for (size_t i = 0; i < a.first.size(); ++i) {
        if (a.first[i] < b.first[i]) return true;
        if (a.first[i] > b.first[i]) return false;
    }

    // If column vectors are identical, compare the size of the value vectors
    if (a.second.size() < b.second.size()) return true;
    if (a.second.size() > b.second.size()) return false;

    // If value vector sizes are the same, compare their elements
    for (size_t i = 0; i < a.second.size(); ++i) {
        if (a.second[i] < b.second[i]) return true;
        if (a.second[i] > b.second[i]) return false;
    }

    // If they are identical, return false
    return false;
}

/**
 * @brief Helper function to get all rows for a d-set (union of its interactions).
 *
 * @param d_set The set of interactions.
 * @param A The covering array.
 * @return A set of row indices that cover *at least one* interaction in 'd_set'.
 */
robin_hood::unordered_flat_set<int> rows_of_d_set(const d_set_type& d_set, const ca_type& A) {
    robin_hood::unordered_flat_set<int> the_rows;
    for(const auto& interaction : d_set) {
        // Get rows for this single interaction
        auto rows = rows_of_interaction(interaction, A);
        // Add them to the union set
        the_rows.insert(rows.begin(), rows.end());
    }
    return the_rows;
}

/**
 * @brief Finds all non-detecting (interaction, d-set) pairs.
 *
 * The detecting condition (Colbourn & McClary 2008, §2):
 *   For every d-set T (|T| ≤ d, independent) and every interaction X ∉ T
 *   (independent of T): |ρ(A,X) \ ρ(A,T)| ≥ λ.
 *
 * A pair (X, T) VIOLATES this when |ρ(A,X) \ ρ(A,T)| < λ, meaning
 * the rows covering X are too heavily contained within the rows of T.
 *
 * OUTPUT: vector of (singleton_dset_id_of_X, dset_id_of_T, separation).
 *   The first element is always a size-1 d-set encoding of interaction X.
 *   The second element is the d-set T.
 *   The third is |ρ(A,X) \ ρ(A,T)| in the initial array.
 *
 * Phase 2 then adds rows where each new row that covers X but does NOT
 * cover any interaction in T contributes +1 to the separation count.
 *
 * Uses partition counts for pruning: the one-sided lower bound
 *   Σ_p max(0, c_X[p] - c_T[p]) ≥ |ρ(A,X) \ ρ(A,T)|
 * allows skipping pairs where this bound ≥ λ.
 *
 * @param codec  [in/out] The codec to initialize and use for encoding.
 * @return A vector of (singleton_dset_id, dset_id, separation) tuples.
 */
std::vector<undist_pair_type> find_non_detecting_sets( const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar, int X, InteractionCodec& codec) {
    
    // 1. Generate all possible interactions and initialize the codec
    auto interactions = get_interactions(t, vs, t_bar);
    codec.init(interactions, d);
    std::cout << "Codec initialized: T=" << codec.T << " interactions, d_max=" << codec.d_max << ", M=" << codec.M << " possible d-sets\n";

    int N = (int)A.size();
    int partition_size = std::max(1, N / X);

    // 2. Build RowBitsets for each interaction (fast set-minus via AND-NOT)
    std::vector<RowBitset> ix_bitsets(codec.T, RowBitset(N));
    for (int r = 0; r < N; ++r)
        for (interaction_id iid = 0; iid < codec.T; ++iid) {
            const auto& ix = codec.decode_interaction(iid);
            bool covers = true;
            for (size_t j = 0; j < ix.first.size(); ++j)
                if (A[r][ix.first[j]] != ix.second[j]) { covers = false; break; }
            if (covers) ix_bitsets[iid].set(r);
        }

    // 3. Build partition counts for each interaction
    std::vector<partitioned_counts> ix_pcounts(codec.T, partitioned_counts(X, 0));
    for (interaction_id iid = 0; iid < codec.T; ++iid) {
        for (int r = 0; r < N; ++r) {
            if (ix_bitsets[iid].test(r)) {
                int pi = std::min(X - 1, r / partition_size);
                ix_pcounts[iid][pi]++;
            }
        }
    }

    // 4. Generate all valid d-sets: their ids, interaction_ids, bitsets, partition counts
    auto lower_lim = d;
    if (d_bar) lower_lim = 1;

    struct DSetInfo {
        d_set_id id;
        std::vector<interaction_id> iids;
        RowBitset bitset;
        partitioned_counts pcounts;
    };
    std::vector<DSetInfo> all_dsets;

    std::cout << "Generating d-sets for detecting array check...\n";

    for (int sz = lower_lim; sz <= d; sz++) {
        auto combos = combinations(interactions, sz);
        for (auto& combo : combos) {
            d_set_type inner_dset;
            for (auto& interaction : combo) inner_dset.push_back(interaction);

            // Subsumption filter
            if (sz > 1 && has_redundant_interaction(inner_dset)) continue;

            DSetInfo info;
            info.id = codec.encode_d_set(inner_dset);

            // Get interaction ids
            info.iids = codec.decode_d_set_ids(info.id);

            // Build bitset (union of interaction bitsets)
            info.bitset = RowBitset(N);
            for (auto iid : info.iids) info.bitset.union_with(ix_bitsets[iid]);

            // Build partition counts from the union bitset
            info.pcounts.assign(X, 0);
            for (int r = 0; r < N; ++r) {
                if (info.bitset.test(r)) {
                    int pi = std::min(X - 1, r / partition_size);
                    info.pcounts[pi]++;
                }
            }

            all_dsets.push_back(std::move(info));
        }
    }

    std::cout << "Generated " << all_dsets.size() << " valid d-sets. "
              << "Checking detecting condition for " << codec.T
              << " interactions...\n";

    // 5. For each (interaction X, d-set T) pair, check the detecting condition
    std::vector<undist_pair_type> to_return;
    size_t pairs_checked = 0;
    size_t pairs_pruned = 0;

    for (interaction_id x_iid = 0; x_iid < codec.T; ++x_iid) {
        // Encode X as a singleton d-set
        d_set_id x_singleton = codec.encode_d_set_ids({x_iid});

        for (size_t ti = 0; ti < all_dsets.size(); ++ti) {
            const auto& dset_info = all_dsets[ti];

            // Skip if X ∈ T
            bool x_in_t = false;
            for (auto iid : dset_info.iids) {
                if (iid == x_iid) { x_in_t = true; break; }
            }
            if (x_in_t) continue;

            // Independence check (for t_bar): X must be independent of T
            // No interaction in T subsumes X, and X doesn't subsume any in T
            if (t_bar) {
                const auto& x_ix = codec.decode_interaction(x_iid);
                bool dependent = false;
                for (auto iid : dset_info.iids) {
                    const auto& t_ix = codec.decode_interaction(iid);
                    if (interaction_subsumes(x_ix, t_ix) ||
                        interaction_subsumes(t_ix, x_ix)) {
                        dependent = true;
                        break;
                    }
                }
                if (dependent) continue;
            }

            // Partition-count pruning: one-sided lower bound
            //   Σ_p max(0, c_X[p] - c_T[p]) is a lower bound on |ρ(X)\ρ(T)|
            int lb = 0;
            for (int p = 0; p < X; ++p) {
                int diff = ix_pcounts[x_iid][p] - dset_info.pcounts[p];
                if (diff > 0) lb += diff;
            }
            if (lb >= lambda) {
                pairs_pruned++;
                continue;  // guaranteed to be distinguished
            }

            // Full check via bitset
            int sep = ix_bitsets[x_iid].andnot_popcount(dset_info.bitset);
            pairs_checked++;

            if (sep < lambda) {
                to_return.push_back({x_singleton, dset_info.id, sep});
            }
        }
    }

    std::cout << "Detecting check complete: " << pairs_checked << " pairs checked, "
              << pairs_pruned << " pruned by partition bound, "
              << to_return.size() << " violations found.\n";

    return to_return;
}

/**
 * @brief Finds all non-locating pairs of d-sets.
 *
 * A pair of d-sets (D1, D2) is non-locating if the symmetric difference
 * of their row sets is less than lambda.
 * | R(D1) XOR R(D2) | < lambda
 *
 * CORRECTNESS FIX: For lambda > 1, exact-match bucketing on partition
 * signatures misses non-locating pairs whose signatures differ slightly.
 *
 * The fix:
 *   1. Compute ACTUAL d-set partition counts using the UNION of rows
 *      (not the sum of interaction counts, which overcounts).
 *   2. The L1 distance between two actual partition count vectors is
 *      a LOWER BOUND on the symmetric difference (since partitions
 *      are disjoint and |A△B| >= ||A|-|B|| per partition).
 *   3. Only pairs with L1 distance < lambda need the expensive full check.
 *
 * @param codec  [in/out] The codec to initialize and use for encoding.
 * @return A vector of (d_set_id, d_set_id, diff_size) tuples.
 */
std::vector<undist_pair_type> find_non_locating_sets(const ca_type& A, t_type t, const vs_type& vs, lambda_type lambda, d_type d, bool d_bar, bool t_bar, int X, InteractionCodec& codec) {
    // 1. Generate all interactions and initialize the codec
    auto interactions = get_interactions(t, vs, t_bar);
    codec.init(interactions, d);
    std::cout << "Codec initialized: T=" << codec.T << " interactions, d_max=" << codec.d_max << ", M=" << codec.M << " possible d-sets\n";

    // 2. Generate all d-sets, encode them, compute ACTUAL partition counts
    auto lower_lim = d;
    if (d_bar) {
        lower_lim = 1;
    }

    // Store each d-set's id and its actual partition count signature
    struct DSetEntry {
        d_set_id id;
        partitioned_counts counts; // actual row counts per partition (via union)
    };
    std::vector<DSetEntry> all_d_sets;

    int partition_size = std::max(1, static_cast<int>(A.size() / X));

    std::cout << "Computing actual partitioned row counts of d-sets...\n";

    for (int i = lower_lim; i <= d; i++) {
        auto to_add = combinations(interactions, i);
        for (auto& individual_d_set : to_add) {
            d_set_type inner_d_set;
            for (auto& interaction : individual_d_set) {
                inner_d_set.push_back(interaction);
            }

            // Skip d-sets with redundant (subsumed) interactions.
            if (i > 1 && has_redundant_interaction(inner_d_set)) continue;

            d_set_id id = codec.encode_d_set(inner_d_set);

            // Compute ACTUAL partition counts using the UNION of rows
            partitioned_counts total_counts(X, 0);
            robin_hood::unordered_set<int> distinct_rows;
            for (const auto& interaction : inner_d_set) {
                auto rows = rows_of_interaction(interaction, A);
                distinct_rows.insert(rows.begin(), rows.end());
            }
            for (const auto& row_idx : distinct_rows) {
                int partition_index = std::min(X - 1, static_cast<int>(row_idx / partition_size));
                total_counts[partition_index]++;
            }

            all_d_sets.push_back({id, std::move(total_counts)});
        }
    }

    // 3. Group by partition signature for fast intra-bucket pairing
    robin_hood::unordered_map<partitioned_counts, std::vector<size_t>, VectorHasher> buckets;
    for (size_t idx = 0; idx < all_d_sets.size(); ++idx) {
        buckets[all_d_sets[idx].counts].push_back(idx);
    }

    // 4. Helper: full symmetric difference check between two d-sets
    auto check_pair = [&](size_t idx_i, size_t idx_j,
                          std::vector<undist_pair_type>& results) {
        d_set_type dset_i = codec.decode_d_set(all_d_sets[idx_i].id);
        d_set_type dset_j = codec.decode_d_set(all_d_sets[idx_j].id);

        auto unordered_rows1 = rows_of_d_set(dset_i, A);
        auto unordered_rows2 = rows_of_d_set(dset_j, A);
        std::vector<int> rows1(unordered_rows1.begin(), unordered_rows1.end());
        std::vector<int> rows2(unordered_rows2.begin(), unordered_rows2.end());
        std::sort(rows1.begin(), rows1.end());
        std::sort(rows2.begin(), rows2.end());

        int diff_size = size_of_symmetric_difference(
            rows1.begin(), rows1.end(), rows2.begin(), rows2.end());

        if (diff_size < lambda) {
            results.push_back(std::make_tuple(
                all_d_sets[idx_i].id, all_d_sets[idx_j].id, diff_size));
        }
    };

    // 5. Pre-enumerate all delta vectors in Z^X with ||delta||_1 <= lambda-1.
    //    These are the L1 ball neighbors.  For each bucket signature s,
    //    we look up s+delta in the hash map — O(B * |ball|) instead of O(B^2).
    int radius = (int)lambda - 1;

    std::vector<partitioned_counts> deltas;
    {
        // Recursive generator: fill delta[pos..X-1] with remaining L1 budget.
        partitioned_counts current(X, 0);
        std::function<void(int, int)> gen_deltas = [&](int pos, int budget) {
            if (pos == X) {
                deltas.push_back(current);
                return;
            }
            // Try delta[pos] = 0, ±1, ±2, ..., ±budget
            for (int v = -budget; v <= budget; ++v) {
                current[pos] = v;
                gen_deltas(pos + 1, budget - std::abs(v));
            }
            current[pos] = 0;
        };
        gen_deltas(0, radius);
    }

    std::cout << "L1 ball size for radius " << radius
              << " in Z^" << X << ": " << deltas.size() << " vectors\n";

    // 6. Find non-locating pairs via neighbor lookup
    std::vector<undist_pair_type> to_return;
    std::cout << "Ready to look at pairs via neighbor enumeration...\n";

    // For deduplication: only process bucket pair (s, s') when s <= s' (lex).
    // When delta is the zero vector, s = s' so we do intra-bucket pairs.
    // Otherwise, we compute s' = s + delta and check s <= s'.
    size_t bucket_pairs_checked = 0;

    for (const auto& [sig, indices] : buckets) {
        for (const auto& delta : deltas) {
            // Compute neighbor signature s' = sig + delta
            partitioned_counts neighbor(X);
            bool valid = true;
            for (int p = 0; p < X; ++p) {
                neighbor[p] = sig[p] + delta[p];
                // Partition counts can't be negative
                if (neighbor[p] < 0) { valid = false; break; }
            }
            if (!valid) continue;

            // Look up in hash map
            auto it = buckets.find(neighbor);
            if (it == buckets.end()) continue;

            // Deduplication: only process when sig <= neighbor (lex)
            if (neighbor < sig) continue;

            const auto& other_indices = it->second;
            bucket_pairs_checked++;

            if (sig == neighbor) {
                // Same bucket (delta = 0): intra-bucket pairs
                for (size_t i = 0; i < indices.size(); ++i) {
                    for (size_t j = i + 1; j < indices.size(); ++j) {
                        check_pair(indices[i], indices[j], to_return);
                    }
                }
            } else {
                // Cross-bucket: all pairs between sig and neighbor
                for (size_t i = 0; i < indices.size(); ++i) {
                    for (size_t j = 0; j < other_indices.size(); ++j) {
                        check_pair(indices[i], other_indices[j], to_return);
                    }
                }
            }
        }
    }

    std::cout << "Checked " << bucket_pairs_checked << " bucket pairs (of "
              << buckets.size() << " total buckets)\n";

    return to_return;
}
