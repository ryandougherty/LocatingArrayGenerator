#pragma once
#include "../LocAG.h"

/**
 * @brief Runs the Conditional Expectation (CE) algorithm from the 2022 paper.
 *
 * This algorithm constructs the array one row at a time. For each row,
 * it greedily selects the value for each factor that minimizes the
 * expected number of interactions that will remain uncovered (i.e., not
 * covered lambda times) by the end of the entire process.
 *
 * @param array A pointer to the LocatingArray object, which contains
 * all parameters (k, v, t, l) and the 
 * 'uncovered_interactions' map. The function will populate
 * 'array->array' with the generated rows.
 */
void run_phase_2_ce(LocatingArray *array);