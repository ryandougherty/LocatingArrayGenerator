#pragma once
#include <random>
#include "../LocAG.h"

struct PercentGAFitnessInd {
    std::vector<double> percents;
    int N = -1;
    long long time = -1;
    ca_type generated_rows;

    bool operator==(PercentGAFitnessInd const&) const = default;
};

std::vector<PercentGAFitnessInd> percent_GA(d_type d, t_type t, const vs_type& vs, const lambda_type& l, const std::vector<undist_pair_type>& non_locating_pairs, bool use_default_percents, bool is_detecting, const std::string& execution_policy, const InteractionCodec& codec);
