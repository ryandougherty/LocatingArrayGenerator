/* ----------------------------------------------------------------------------
 * phase2.cpp
 *
 * Phase 2: optimization stage (GA / SA).
 *
 * OPTIMIZATIONS APPLIED:
 *   1. RowBitset — row sets stored as bitsets; symmetric difference is
 *      XOR + popcount instead of O(N) merge scan.
 *   2. Precomputed interaction→RowBitset table — built once per candidate
 *      array, d-set bitsets derived via OR.
 *   3. Incremental SA — on each mutation, only the affected interactions,
 *      d-sets, and pairs are recomputed.  Accept/reject via undo.
 * ----------------------------------------------------------------------------
 */

#include "../utils/utils.h"
#include "../utils/row_bitset.h"
#include "phase2.h"
#include "../phase1/phase1.h"

#ifdef HAS_TBB
#include <execution>
#endif

/* ============================= helpers =================================== */

/** Does row r of array A cover interaction ix? */
static inline bool row_covers_ix(const std::vector<v_type>& row,
                                 const interaction_type& ix) {
    for (size_t j = 0; j < ix.first.size(); ++j)
        if (row[ix.first[j]] != ix.second[j]) return false;
    return true;
}

/** Build interaction→RowBitset table for an entire array. O(N·T·t) */
static std::vector<RowBitset> build_interaction_bitsets(
        const ca_type& A, int N, const InteractionCodec& codec) {
    std::vector<RowBitset> irows(codec.T, RowBitset(N));
    for (int r = 0; r < N; ++r)
        for (interaction_id iid = 0; iid < codec.T; ++iid)
            if (row_covers_ix(A[r], codec.decode_interaction(iid)))
                irows[iid].set(r);
    return irows;
}

/* ==================== static index structures ============================ */

/**
 * SAIndexes: built ONCE per set of pairs, reused across every SA call
 * in go().  Contains the column→interaction, interaction→d-set, and
 * d-set→pair index mappings.
 */
struct SAIndexes {
    /* col → interaction_ids involving that column */
    std::vector<std::vector<interaction_id>> col_to_iids;

    /* For each unique d-set in the pair list: its constituent iids */
    std::unordered_map<d_set_id, std::vector<interaction_id>> dset_iids;

    /* iid → d-set ids containing it (restricted to d-sets in pairs) */
    std::vector<std::vector<d_set_id>> iid_to_dsets;

    /* d-set id → indices into the pair vector */
    std::unordered_map<d_set_id, std::vector<size_t>> dset_to_pairs;

    /* Per-pair requirement: how much MORE separation is needed */
    std::vector<int> requirements;

    void build(int k,
               const InteractionCodec& codec,
               const std::vector<undist_pair_type>& pairs,
               bool is_detecting, lambda_type l)
    {
        /* 1. col → interactions */
        col_to_iids.assign(k, {});
        for (interaction_id iid = 0; iid < codec.T; ++iid)
            for (auto c : codec.decode_interaction(iid).first)
                col_to_iids[c].push_back(iid);

        /* 2. unique d-sets → constituent iids */
        dset_iids.clear();
        for (auto& [id1, id2, sep] : pairs) {
            if (!dset_iids.count(id1))
                dset_iids[id1] = codec.decode_d_set_ids(id1);
            if (!dset_iids.count(id2))
                dset_iids[id2] = codec.decode_d_set_ids(id2);
        }

        /* 3. iid → d-sets */
        iid_to_dsets.assign(codec.T, {});
        for (auto& [did, iids] : dset_iids)
            for (auto iid : iids)
                iid_to_dsets[iid].push_back(did);

        /* 4. d-set → pair indices */
        dset_to_pairs.clear();
        for (size_t pi = 0; pi < pairs.size(); ++pi) {
            auto& [id1, id2, sep] = pairs[pi];
            dset_to_pairs[id1].push_back(pi);
            dset_to_pairs[id2].push_back(pi);
        }

        /* 5. per-pair requirement */
        requirements.resize(pairs.size());
        for (size_t pi = 0; pi < pairs.size(); ++pi) {
            auto& [id1, id2, sep] = pairs[pi];
            requirements[pi] = (int)l - sep;
        }
    }
};

/* ========================= mutable SA state ============================== */

/**
 * SAState: holds the current candidate array and all derived bitset
 * structures.  Supports apply_mutation / undo_mutation for the
 * incremental SA.
 */
struct SAState {
    ca_type A;
    int N;
    bool is_detecting = false;

    std::vector<RowBitset>                     irows;      /* T entries */
    std::unordered_map<d_set_id, RowBitset>    drows;
    std::vector<int>                           pair_seps;  /* per-pair symm-diff */
    int                                        score;

    /* ---------- initialise from a fresh array ----------------------------- */
    void init(const ca_type& arr,
              const SAIndexes& idx,
              const InteractionCodec& codec,
              const std::vector<undist_pair_type>& pairs,
              bool detecting)
    {
        A = arr;
        N = (int)A.size();
        is_detecting = detecting;

        /* interaction bitsets */
        irows = build_interaction_bitsets(A, N, codec);

        /* d-set bitsets (OR of constituent interaction bitsets) */
        drows.clear();
        for (auto& [did, iids] : idx.dset_iids) {
            RowBitset bs(N);
            for (auto iid : iids) bs.union_with(irows[iid]);
            drows[did] = std::move(bs);
        }

        /* pair separations & score */
        pair_seps.resize(pairs.size());
        score = 0;
        for (size_t pi = 0; pi < pairs.size(); ++pi) {
            auto& [id1, id2, sep] = pairs[pi];
            pair_seps[pi] = drows.at(id1).separation(drows.at(id2), is_detecting);
            if (pair_seps[pi] >= idx.requirements[pi]) ++score;
        }
    }

    /* ---------- mutation + undo machinery -------------------------------- */

    /* saved state for undo */
    struct Undo {
        int mut_row, mut_col;
        std::vector<v_type> old_vals;                              /* array cells */
        std::vector<std::pair<interaction_id, RowBitset>> old_irows;
        std::vector<std::pair<d_set_id, RowBitset>>      old_drows;
        std::vector<std::pair<size_t, int>>               old_pair_seps;
        int old_score;
    };
    Undo undo;

    /**
     * Apply a random mutation, track changes, return the new score.
     *   mut_type: 0 = row, 1 = column, 2 = cell  (matching original mutate())
     */
    int apply_mutation(int mut_type, std::mt19937& rng_local,
                       const SAIndexes& idx, const InteractionCodec& codec,
                       const std::vector<undist_pair_type>& pairs,
                       const vs_type& vs)
    {
        undo.old_score = score;
        undo.old_irows.clear();
        undo.old_drows.clear();
        undo.old_pair_seps.clear();

        /* ---- choose what to mutate ---- */

        if (mut_type == 0) {
            /* mutate entire row */
            int r = any_int(rng_local) % N;
            undo.mut_row = r; undo.mut_col = -1;
            undo.old_vals.assign(A[r].begin(), A[r].end());
            for (size_t c = 0; c < vs.size(); ++c)
                A[r][c] = any_int(rng_local) % vs[c];
            update_after_row_change(r, idx, codec, pairs);

        } else if (mut_type == 1) {
            /* mutate entire column */
            int c = any_int(rng_local) % (int)vs.size();
            undo.mut_row = -1; undo.mut_col = c;
            undo.old_vals.resize(N);
            for (int r = 0; r < N; ++r) {
                undo.old_vals[r] = A[r][c];
                A[r][c] = any_int(rng_local) % vs[c];
            }
            update_after_col_change(c, idx, codec, pairs);

        } else {
            /* mutate single cell */
            int r = any_int(rng_local) % N;
            int c = any_int(rng_local) % (int)vs.size();
            undo.mut_row = r; undo.mut_col = c;
            undo.old_vals = { A[r][c] };
            A[r][c] = any_int(rng_local) % vs[c];
            update_after_cell_change(r, c, idx, codec, pairs);
        }

        return score;
    }

    void reject_mutation() {
        /* restore array */
        if (undo.mut_col == -1 && undo.mut_row >= 0) {
            /* row mutation */
            A[undo.mut_row] = undo.old_vals;
        } else if (undo.mut_row == -1 && undo.mut_col >= 0) {
            /* col mutation */
            for (int r = 0; r < N; ++r)
                A[r][undo.mut_col] = undo.old_vals[r];
        } else {
            /* cell mutation */
            A[undo.mut_row][undo.mut_col] = undo.old_vals[0];
        }
        /* restore bitsets */
        for (auto& [iid, old_bs] : undo.old_irows) irows[iid] = std::move(old_bs);
        for (auto& [did, old_bs] : undo.old_drows) drows[did] = std::move(old_bs);
        for (auto& [pi, old_sep] : undo.old_pair_seps) pair_seps[pi] = old_sep;
        score = undo.old_score;
    }

    /* accept = do nothing (changes already in place) */
    void accept_mutation() { }

private:

    /* ---- propagation after interaction bitset changes ---- */

    void propagate(const std::vector<interaction_id>& changed_iids,
                   const SAIndexes& idx,
                   const std::vector<undist_pair_type>& pairs)
    {
        /* collect affected d-sets */
        std::unordered_set<d_set_id> affected_dsets;
        for (auto iid : changed_iids)
            for (auto did : idx.iid_to_dsets[iid])
                affected_dsets.insert(did);

        /* recompute affected d-set bitsets */
        for (auto did : affected_dsets) {
            undo.old_drows.push_back({did, drows[did]});
            RowBitset bs(N);
            for (auto iid : idx.dset_iids.at(did))
                bs.union_with(irows[iid]);
            drows[did] = std::move(bs);
        }

        /* collect affected pairs */
        std::unordered_set<size_t> affected_pairs;
        for (auto did : affected_dsets)
            if (idx.dset_to_pairs.count(did))
                for (auto pi : idx.dset_to_pairs.at(did))
                    affected_pairs.insert(pi);

        /* recompute affected pair separations & update score */
        for (auto pi : affected_pairs) {
            undo.old_pair_seps.push_back({pi, pair_seps[pi]});
            bool was_ok = (pair_seps[pi] >= idx.requirements[pi]);
            auto& [id1, id2, sep] = pairs[pi];
            pair_seps[pi] = drows.at(id1).separation(drows.at(id2), is_detecting);
            bool now_ok = (pair_seps[pi] >= idx.requirements[pi]);
            if (was_ok && !now_ok) --score;
            else if (!was_ok && now_ok) ++score;
        }
    }

    /* ---- cell change: only interactions involving column c, row r ---- */

    void update_after_cell_change(int r, int c,
            const SAIndexes& idx, const InteractionCodec& codec,
            const std::vector<undist_pair_type>& pairs)
    {
        std::vector<interaction_id> changed;
        for (auto iid : idx.col_to_iids[c]) {
            bool was = irows[iid].test(r);
            bool now = row_covers_ix(A[r], codec.decode_interaction(iid));
            if (was != now) {
                undo.old_irows.push_back({iid, irows[iid]});
                if (now) irows[iid].set(r); else irows[iid].clear(r);
                changed.push_back(iid);
            }
        }
        if (!changed.empty()) propagate(changed, idx, pairs);
    }

    /* ---- row change: check all T interactions for row r ---- */

    void update_after_row_change(int r,
            const SAIndexes& idx, const InteractionCodec& codec,
            const std::vector<undist_pair_type>& pairs)
    {
        std::vector<interaction_id> changed;
        for (interaction_id iid = 0; iid < (interaction_id)codec.T; ++iid) {
            bool was = irows[iid].test(r);
            bool now = row_covers_ix(A[r], codec.decode_interaction(iid));
            if (was != now) {
                undo.old_irows.push_back({iid, irows[iid]});
                if (now) irows[iid].set(r); else irows[iid].clear(r);
                changed.push_back(iid);
            }
        }
        if (!changed.empty()) propagate(changed, idx, pairs);
    }

    /* ---- column change: rebuild affected interaction bitsets ---- */

    void update_after_col_change(int c,
            const SAIndexes& idx, const InteractionCodec& codec,
            const std::vector<undist_pair_type>& pairs)
    {
        std::vector<interaction_id> changed;
        for (auto iid : idx.col_to_iids[c]) {
            undo.old_irows.push_back({iid, irows[iid]});
            /* rebuild this interaction's bitset from scratch */
            RowBitset bs(N);
            const auto& ix = codec.decode_interaction(iid);
            for (int r = 0; r < N; ++r)
                if (row_covers_ix(A[r], ix))
                    bs.set(r);
            if (bs != irows[iid]) {
                irows[iid] = std::move(bs);
                changed.push_back(iid);
            }
        }
        if (!changed.empty()) propagate(changed, idx, pairs);
    }
};

/* ========================== bitset-based fitness ========================= */

/**
 * Bitset-based fitness — used by the (deprecated) GA path.
 * Builds the full bitset tables from scratch for a candidate array.
 */
int fitness(const ca_type& ind, d_type d, t_type t, const vs_type& vs,
            lambda_type l,
            const std::vector<undist_pair_type>& non_locating_pairs,
            const int& threshold, bool is_detecting,
            const InteractionCodec& codec)
{
    int N = (int)ind.size();
    auto irows = build_interaction_bitsets(ind, N, codec);

    /* d-set bitset cache (integer keys, no custom hasher) */
    std::unordered_map<d_set_id, RowBitset> drows;
    auto get_drow = [&](d_set_id did) -> const RowBitset& {
        auto it = drows.find(did);
        if (it != drows.end()) return it->second;
        auto iids = codec.decode_d_set_ids(did);
        RowBitset bs(N);
        for (auto iid : iids) bs.union_with(irows[iid]);
        drows[did] = std::move(bs);
        return drows[did];
    };

    int score = 0;
    for (auto& [id1, id2, already] : non_locating_pairs) {
        int req = l - already;
        int n = get_drow(id1).separation(get_drow(id2), is_detecting);
        if (n >= req) ++score;
        if (score >= threshold) return threshold + 1;
    }
    return score;
}

/* ========================= GA operators (unchanged) ====================== */

ca_type cross(const ca_type& p1, const ca_type& p2,
              d_type d, t_type t, const vs_type& vs, lambda_type l,
              std::mt19937& rng_local) {
    int val = any_int(rng_local) % 2;
    int n = (int)p1.size();
    ca_type child;
    if (val == 0) {
        auto ri = any_int(rng_local) % p1.size();
        for (size_t i = 0; i < ri; i++) child.push_back(p1[i]);
        for (size_t i = ri; i < (size_t)n; i++) child.push_back(p2[i]);
    } else if (p1.size() != 1) {
        auto r1 = any_int(rng_local) % p1.size();
        auto r2 = any_int(rng_local) % p1.size();
        while (r1 == r2) r2 = any_int(rng_local) % p1.size();
        auto lo = std::min(r1, r2), hi = std::max(r1, r2);
        for (size_t i = 0; i < lo; i++)  child.push_back(p1[i]);
        for (size_t i = lo; i < hi; i++) child.push_back(p2[i]);
        for (size_t i = hi; i < p1.size(); i++) child.push_back(p1[i]);
    } else {
        child = (any_int(rng_local) % 2 == 0) ? p1 : p2;
    }
    return child;
}

ca_type mutate(const ca_type& p1, d_type d, t_type t, const vs_type& vs,
               lambda_type l, std::mt19937& rng_local) {
    int val = any_int(rng_local) % 3;
    int n = (int)p1.size();
    ca_type child = p1;
    if (val == 0) {
        int r = any_int(rng_local) % n;
        for (size_t c = 0; c < vs.size(); c++)
            child[r][c] = any_int(rng_local) % vs[c];
    } else if (val == 1) {
        int c = any_int(rng_local) % (int)vs.size();
        for (int r = 0; r < n; r++)
            child[r][c] = any_int(rng_local) % vs[c];
    } else {
        int r = any_int(rng_local) % n;
        int c = any_int(rng_local) % (int)vs.size();
        child[r][c] = any_int(rng_local) % vs[c];
    }
    return child;
}

struct Ind_NonRecompute_Fitness { ca_type A; int fitness; };

/* ==================== try_N (deprecated GA, uses bitset fitness) ========= */

ca_type try_N(N_type N, d_type d, t_type t, const vs_type& vs, lambda_type l,
              const std::vector<undist_pair_type>& pairs, double percent,
              bool is_detecting, std::mt19937& rng_local,
              const InteractionCodec& codec) {
    ca_type s;
    int pop_size = 100, num_gens = 50;
    std::vector<Ind_NonRecompute_Fitness> pop(pop_size);
    for (auto& e : pop) { e.A = random_array(N, vs.size(), vs); e.fitness = -1; }
    const int max_f = (int)(pairs.size() * percent);
    for (int gen = 0; gen < num_gens; gen++) {
        std::vector<std::pair<int, Ind_NonRecompute_Fitness>> fits;
        for (auto& I : pop) {
            int f = I.fitness;
            if (f == -1) { f = fitness(I.A,d,t,vs,l,pairs,max_f,is_detecting,codec); I.fitness = f; }
            if (f >= max_f) return I.A;
            fits.push_back({f, I});
        }
        std::sort(fits.begin(), fits.end(), [](auto& a, auto& b){ return a.first < b.first; });
        std::vector<Ind_NonRecompute_Fitness> nv;
        for (int i = pop_size/2; i < pop_size; i++) nv.push_back(fits[i].second);
        pop = nv; nv.clear();
        while (nv.size() < (size_t)pop_size/2) {
            auto& p1 = pop[any_int(rng_local) % pop.size()];
            auto& p2 = pop[any_int(rng_local) % pop.size()];
            if (any_int(rng_local) % 10 == 0) {
                auto c = cross(p1.A, p2.A, d,t,vs,l,rng_local);
                if (any_int(rng_local) % 10 < 3) c = mutate(c,d,t,vs,l,rng_local);
                nv.push_back({c, -1});
            }
        }
        pop.insert(pop.end(), nv.begin(), nv.end());
    }
    return s;
}

/* ============= try_N_SA — INCREMENTAL Simulated Annealing ================ */

/**
 * Tries to find an array of N rows that fixes all pairs.
 * Uses the SAState + SAIndexes infrastructure for O(affected) updates
 * instead of O(pairs·d·N) per iteration.
 */
ca_type try_N_SA(N_type N, d_type d, t_type t, const vs_type& vs,
                 lambda_type l,
                 const std::vector<undist_pair_type>& pairs,
                 bool is_detecting, std::mt19937& rng_local,
                 const InteractionCodec& codec,
                 const SAIndexes& sa_idx)
{
    ca_type empty;

    /* initialise SA state */
    SAState state;
    state.init(random_array(N, vs.size(), vs), sa_idx, codec, pairs, is_detecting);

    int required = (int)pairs.size();
    if (state.score >= required) return state.A;

    double temp = 1.0, rate = 0.99;
    int num_iter = 1000;
    int f = state.score;

    for (int it = 0; it < num_iter; ++it) {
        if (f >= required) return state.A;

        int mut_type = any_int(rng_local) % 3;
        int f_new = state.apply_mutation(mut_type, rng_local, sa_idx, codec, pairs, vs);

        if (f_new >= f) {
            state.accept_mutation();
            f = f_new;
        } else {
            double prob = std::exp((double)(f_new - f) / temp);
            if (prob > unif(rng_local)) {
                state.accept_mutation();
                f = f_new;
            } else {
                state.reject_mutation();
                /* f unchanged */
            }
        }
        temp *= rate;
    }

    return (f >= required) ? state.A : empty;
}

/* ========================= go() — binary-search wrapper ================== */

/**
 * Finds the minimal N to fix a percentage of pairs.
 * Builds SAIndexes ONCE and reuses across all SA calls.
 */
ca_type go(const d_type& d, const t_type& t, const vs_type& vs,
           const lambda_type& l,
           const std::vector<undist_pair_type>& non_locating_pairs,
           const double& percent, bool is_detecting,
           std::mt19937& rng_local,
           const InteractionCodec& codec)
{
    bool succ_first = true;
    ca_type result;

    const std::vector<undist_pair_type> only(
        non_locating_pairs.begin(),
        non_locating_pairs.begin() + (size_t)(non_locating_pairs.size() * percent));

    /* build static SA indexes once for this pair subset */
    SAIndexes sa_idx;
    sa_idx.build((int)vs.size(), codec, only, is_detecting, l);

    int N = 1;
    for (auto& [d1, d2, num] : only) N = std::max(N, (int)(l - num));

    /* exponential search */
    while (true) {
        result = try_N_SA(N, d, t, vs, l, only, is_detecting, rng_local, codec, sa_idx);
        if (succ_first && result.size() > 0) return result;
        if (result.size() > 0) break;
        N *= 2;
        succ_first = false;
    }

    /* binary search */
    int N_hi = N, N_lo = N / 2;
    while (N_lo < N_hi) {
        int N_mid = (N_lo + N_hi) / 2;
        auto r2 = try_N_SA(N_mid, d, t, vs, l, only, is_detecting, rng_local, codec, sa_idx);
        if (r2.size() > 0) { N_hi = N_mid; result = r2; }
        else                  N_lo = N_mid + 1;
    }
    return result;
}

/* ========================= Pareto helpers (unchanged) ==================== */

bool dominates(const PercentGAFitnessInd& a, const PercentGAFitnessInd& b) {
    return a.N <= b.N && a.time <= b.time;
}

auto pareto_and_rest(std::vector<PercentGAFitnessInd> pts) {
    std::vector<PercentGAFitnessInd> dom, par;
    int ci = 0;
    while (true) {
        auto cand = pts[ci]; pts.erase(pts.begin() + ci);
        bool nd = true; size_t i = 0;
        while (!pts.empty() && i < pts.size()) {
            if (dominates(cand, pts[i]))      { dom.push_back(pts[i]); pts.erase(pts.begin()+i); }
            else if (dominates(pts[i], cand)) { nd = false; dom.push_back(cand); i++; }
            else                                i++;
        }
        if (nd) par.push_back(cand);
        if (pts.empty()) break;
    }
    return std::make_pair(par, dom);
}

auto generate_rand_percent_individual() {
    std::vector<double> p;
    for (int i = 0, n = ind_size(rng); i < n; i++) p.push_back(unif(rng));
    std::sort(p.begin(), p.end());
    p[0] = 0.001; p.push_back(1.0);
    PercentGAFitnessInd r; r.percents = p; return r;
}

/* ==================== helper: compute rows for update loop =============== */

static std::vector<N_type> compute_rows_of_dset(
        d_set_id id, const InteractionCodec& codec, const ca_type& arr) {
    d_set_type dset = codec.decode_d_set(id);
    robin_hood::unordered_set<N_type> rs;
    for (auto& ix : dset) {
        auto rows = rows_of_interaction(ix, arr);
        rs.insert(rows.begin(), rows.end());
    }
    std::vector<N_type> v(rs.begin(), rs.end());
    std::sort(v.begin(), v.end());
    return v;
}

/* ==================== run_ga_with_policy ================================= */

template<typename Policy>
auto run_ga_with_policy(
    Policy policy, d_type d, t_type t, const vs_type& vs, lambda_type l,
    const std::vector<undist_pair_type>& non_locating_pairs,
    bool use_default_percents, bool is_detecting,
    const InteractionCodec& codec)
{
    if (use_default_percents) {
        std::mt19937 main_rng(std::random_device{}());
        const std::vector<double> percents = {0.021576,0.021576,0.022644,0.030792,
            0.090424,0.071014,0.083679,0.172455,0.220123,0.415283,1.0};
        auto pairs_copy = non_locating_pairs;
        int num_rows = 0;
        ca_type all_rows;
        auto start = high_resolution_clock::now();
        for (auto pct : percents) {
            std::vector<undist_pair_type> next;
            auto ga = go(d,t,vs,l,pairs_copy,pct,is_detecting,main_rng,codec);
            num_rows += (int)ga.size();
            all_rows.insert(all_rows.end(), ga.begin(), ga.end());
            std::unordered_map<d_set_id, std::vector<N_type>> cache;
            auto get = [&](d_set_id id) -> std::vector<N_type> {
                auto it = cache.find(id); if (it != cache.end()) return it->second;
                auto v = compute_rows_of_dset(id, codec, ga); cache[id]=v; return v;
            };
            for (auto& [d1,d2,sep] : pairs_copy) {
                auto r1=get(d1), r2=get(d2);
                int n = is_detecting
                    ? size_of_set_minus(r1.begin(),r1.end(),r2.begin(),r2.end())
                    : size_of_symmetric_difference(r1.begin(),r1.end(),r2.begin(),r2.end());
                int req = l;
                if (sep + n < req) next.push_back({d1,d2,sep+n});
            }
            pairs_copy = next;
            if (pairs_copy.empty()) break;
            std::cout << "Added " << num_rows << " rows, " << pairs_copy.size() << " pairs remaining\n";
        }
        auto stop = high_resolution_clock::now();
        PercentGAFitnessInd ind;
        ind.N = num_rows; ind.percents = percents;
        ind.time = duration_cast<milliseconds>(stop-start).count();
        ind.generated_rows = all_rows;
        return std::vector<PercentGAFitnessInd>{ind};
    }

    /* ---- full GA ---- */
    int pop_size = 100, num_gens = 50;
    std::vector<PercentGAFitnessInd> result, pop;
    for (int i = 0; i < pop_size; i++) pop.push_back(generate_rand_percent_individual());

    for (int gen = 0; gen < num_gens; gen++) {
        std::cout << "Generation #" << gen << "\n";
        std::for_each(
#ifdef HAS_TBB
            policy,
#endif
            pop.begin(), pop.end(), [&](PercentGAFitnessInd& I) {
            if (I.N != -1 && I.time != -1) return;
            std::mt19937 thr(std::random_device{}());
            long long nr = 0; ca_type all;
            auto pc = non_locating_pairs;
            auto t0 = high_resolution_clock::now();
            for (auto pct : I.percents) {
                std::vector<undist_pair_type> next;
                auto ga = go(d,t,vs,l,pc,pct,is_detecting,thr,codec);
                nr += (long long)ga.size();
                all.insert(all.end(), ga.begin(), ga.end());
                std::unordered_map<d_set_id, std::vector<N_type>> cache;
                auto get = [&](d_set_id id) -> std::vector<N_type> {
                    auto it = cache.find(id); if (it != cache.end()) return it->second;
                    auto v = compute_rows_of_dset(id, codec, ga); cache[id]=v; return v;
                };
                for (auto& [d1,d2,sep] : pc) {
                    auto r1=get(d1), r2=get(d2);
                    int n = is_detecting
                        ? size_of_set_minus(r1.begin(),r1.end(),r2.begin(),r2.end())
                        : size_of_symmetric_difference(r1.begin(),r1.end(),r2.begin(),r2.end());
                    int req = l;
                    if (sep + n < req) next.push_back({d1,d2,sep+n});
                }
                pc = next;
            }
            auto t1 = high_resolution_clock::now();
            I.N = (int)nr; I.time = duration_cast<milliseconds>(t1-t0).count();
            I.generated_rows = all;
        });

        auto fits = pop;
        auto target = (size_t)pop_size/2;
        auto [par, rest] = pareto_and_rest(fits);
        for (auto& x : par) { std::cout << x.N << "," << x.time << ","; print_vec(x.percents); std::cout << "\n"; }
        std::vector<PercentGAFitnessInd> np(par.begin(), par.end());
        result = par;
        while (np.size() < target) {
            for (auto& e : par) fits.erase(std::remove(fits.begin(),fits.end(),e),fits.end());
            auto [p2, r2] = pareto_and_rest(fits);
            for (auto& e : p2) { if (np.size() < target) np.push_back(e); else break; }
            par = p2;
        }
        pop.clear(); for (auto& e : np) pop.push_back(e);

        std::vector<PercentGAFitnessInd> children;
        while (children.size() < (size_t)pop_size/2) {
            auto& p1 = pop[any_int(rng)%pop.size()], &p2 = pop[any_int(rng)%pop.size()];
            auto ri = any_int(rng) % std::min(p1.percents.size(),p2.percents.size());
            while (ri == 0) ri = any_int(rng) % std::min(p1.percents.size(),p2.percents.size());
            PercentGAFitnessInd ch;
            for (size_t i=0; i<ri; i++) ch.percents.push_back(p1.percents[i]);
            for (size_t i=ri; i<p2.percents.size(); i++) ch.percents.push_back(p2.percents[i]);
            children.push_back(ch);
        }
        for (auto& ch : children) {
            int r = any_int(rng)%10;
            if (r==0) { auto i=any_int(rng)%(ch.percents.size()-1); ch.percents.insert(ch.percents.begin()+i, unif(rng)); }
            else if (r==1) { auto i=any_int(rng)%(ch.percents.size()-1); ch.percents.erase(ch.percents.begin()+i); }
            else if (r==2) { auto i=any_int(rng)%(ch.percents.size()-1); ch.percents[i]=unif(rng); }
            pop.push_back(ch);
        }
    }
    return result;
}

/* =========================== public API ================================== */

std::vector<PercentGAFitnessInd> percent_GA(
    d_type d, t_type t, const vs_type& vs, const lambda_type& l,
    const std::vector<undist_pair_type>& non_locating_pairs,
    bool use_default_percents, bool is_detecting,
    const std::string& execution_policy,
    const InteractionCodec& codec)
{
#ifdef HAS_TBB
    if (execution_policy == "parallel") {
        std::cout << "Running GA with std::execution::par\n";
        return run_ga_with_policy(std::execution::par, d,t,vs,l,non_locating_pairs,
                                  use_default_percents, is_detecting, codec);
    } else {
        std::cout << "Running GA with std::execution::seq\n";
        return run_ga_with_policy(std::execution::seq, d,t,vs,l,non_locating_pairs,
                                  use_default_percents, is_detecting, codec);
    }
#else
    if (execution_policy == "parallel")
        std::cout << "WARNING: Parallel not available (no TBB). Falling back to serial.\n";
    std::cout << "Running GA (serial)\n";
    return run_ga_with_policy(0, d,t,vs,l,non_locating_pairs,
                              use_default_percents, is_detecting, codec);
#endif
}
