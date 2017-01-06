//////////////////////////////////////////////////////////////////////
// Config.hpp
//
// Global configuration file.
//////////////////////////////////////////////////////////////////////

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

//////////////////////////////////////////////////////////////////////
// Options related to general inference
//////////////////////////////////////////////////////////////////////

// showing timings for inference routines
#define SHOW_TIMINGS                               0

// use straightforward calculation for FM2 matrix
#define SIMPLE_FM2                                 0

// use candidate list optimization for Viterbi parsing
#define CANDIDATE_LIST                             1

// use unrolled computation for single branch loops
#define FAST_SINGLE_BRANCH_LOOPS                   1

// use caching algorithm for fast helix length scores
#define FAST_HELIX_LENGTHS                         1

//////////////////////////////////////////////////////////////////////
// (A) Options related to max-margin training
//////////////////////////////////////////////////////////////////////

// the maximum loss DELTA(y,y') allocated to each training example; if
// this symbol is undefined, then the DELTA(y,y') loss function is not
// included
//    -- for a straight CRF, this value should be undefined
//    -- for a max-margin model, this value should be set to 1
#define HAMMING_LOSS                               1

//////////////////////////////////////////////////////////////////////
// (E) Used parameter groups
//////////////////////////////////////////////////////////////////////

#ifdef CONTRAFOLD_MODEL

#define PARAMS_BASE_PAIR                           1
#define PARAMS_BASE_PAIR_DIST                      0
#define PARAMS_TERMINAL_MISMATCH                   1
#define PARAMS_HAIRPIN_LENGTH                      1
#define PARAMS_HAIRPIN_3_NUCLEOTIDES               0
#define PARAMS_HAIRPIN_4_NUCLEOTIDES               0
#define PARAMS_HELIX_LENGTH                        0
#define PARAMS_ISOLATED_BASE_PAIR                  0
#define PARAMS_INTERNAL_EXPLICIT                   1
#define PARAMS_BULGE_LENGTH                        1
#define PARAMS_INTERNAL_LENGTH                     1
#define PARAMS_INTERNAL_SYMMETRY                   1
#define PARAMS_INTERNAL_ASYMMETRY                  1
#define PARAMS_BULGE_0x1_NUCLEOTIDES               1
#define PARAMS_BULGE_0x2_NUCLEOTIDES               0
#define PARAMS_BULGE_0x3_NUCLEOTIDES               0
#define PARAMS_INTERNAL_1x1_NUCLEOTIDES            1
#define PARAMS_INTERNAL_1x2_NUCLEOTIDES            0
#define PARAMS_INTERNAL_2x2_NUCLEOTIDES            0
#define PARAMS_HELIX_STACKING                      1
#define PARAMS_HELIX_CLOSING                       1
#define PARAMS_MULTI_LENGTH                        1
#define PARAMS_DANGLE                              1
#define PARAMS_EXTERNAL_LENGTH                     1

#else  // CONTRAFOLD_MODEL

#define PARAMS_BASE_PAIR                           1
#define PARAMS_BASE_PAIR_DIST                      1
#define PARAMS_TERMINAL_MISMATCH                   1
#define PARAMS_HAIRPIN_LENGTH                      1
#define PARAMS_HAIRPIN_3_NUCLEOTIDES               1
#define PARAMS_HAIRPIN_4_NUCLEOTIDES               1
#define PARAMS_HELIX_LENGTH                        1
#define PARAMS_ISOLATED_BASE_PAIR                  1
#define PARAMS_INTERNAL_EXPLICIT                   1
#define PARAMS_BULGE_LENGTH                        1
#define PARAMS_INTERNAL_LENGTH                     1
#define PARAMS_INTERNAL_SYMMETRY                   1
#define PARAMS_INTERNAL_ASYMMETRY                  1
#define PARAMS_BULGE_0x1_NUCLEOTIDES               1
#define PARAMS_BULGE_0x2_NUCLEOTIDES               1
#define PARAMS_BULGE_0x3_NUCLEOTIDES               1
#define PARAMS_INTERNAL_1x1_NUCLEOTIDES            1
#define PARAMS_INTERNAL_1x2_NUCLEOTIDES            1
#define PARAMS_INTERNAL_2x2_NUCLEOTIDES            1
#define PARAMS_HELIX_STACKING                      1
#define PARAMS_HELIX_CLOSING                       1
#define PARAMS_MULTI_LENGTH                        1
#define PARAMS_DANGLE                              1
#define PARAMS_EXTERNAL_LENGTH                     1

// extended parameters
#define PARAMS_HAIRPIN_5_NUCLEOTIDES               1
#define PARAMS_HAIRPIN_6_NUCLEOTIDES               1
#define PARAMS_HAIRPIN_7_NUCLEOTIDES               0

#define PARAMS_BULGE_0x4_NUCLEOTIDES               1
#define PARAMS_BULGE_0x5_NUCLEOTIDES               1
#define PARAMS_BULGE_0x6_NUCLEOTIDES               1
#define PARAMS_BULGE_0x7_NUCLEOTIDES               0

#define PARAMS_INTERNAL_1x3_NUCLEOTIDES            1
#define PARAMS_INTERNAL_2x3_NUCLEOTIDES            1
#define PARAMS_INTERNAL_3x3_NUCLEOTIDES            1

#define PARAMS_INTERNAL_1x4_NUCLEOTIDES            1
#define PARAMS_INTERNAL_2x4_NUCLEOTIDES            1
#define PARAMS_INTERNAL_3x4_NUCLEOTIDES            1
#define PARAMS_INTERNAL_4x4_NUCLEOTIDES            1

#endif  // CONTRAFOLD_MODEL

/*
#define PARAMS_BASE_PAIR                           1
#define PARAMS_BASE_PAIR_DIST                      0
#define PARAMS_TERMINAL_MISMATCH                   0
#define PARAMS_HAIRPIN_LENGTH                      0
#define PARAMS_HAIRPIN_3_NUCLEOTIDES               0
#define PARAMS_HAIRPIN_4_NUCLEOTIDES               0
#define PARAMS_HELIX_LENGTH                        0
#define PARAMS_ISOLATED_BASE_PAIR                  0
#define PARAMS_INTERNAL_EXPLICIT                   0
#define PARAMS_BULGE_LENGTH                        0
#define PARAMS_INTERNAL_LENGTH                     0
#define PARAMS_INTERNAL_SYMMETRY                   0
#define PARAMS_INTERNAL_ASYMMETRY                  0
#define PARAMS_BULGE_0x1_NUCLEOTIDES               0
#define PARAMS_BULGE_0x2_NUCLEOTIDES               0
#define PARAMS_BULGE_0x3_NUCLEOTIDES               0
#define PARAMS_INTERNAL_1x1_NUCLEOTIDES            0
#define PARAMS_INTERNAL_1x2_NUCLEOTIDES            0
#define PARAMS_INTERNAL_2x2_NUCLEOTIDES            0
#define PARAMS_HELIX_STACKING                      0
#define PARAMS_HELIX_CLOSING                       0
#define PARAMS_MULTI_LENGTH                        0
#define PARAMS_DANGLE                              0
#define PARAMS_EXTERNAL_LENGTH                     0
*/

//////////////////////////////////////////////////////////////////////
// (F) Miscellaneous model constants
//////////////////////////////////////////////////////////////////////

const int C_MIN_HAIRPIN_LENGTH = 0;
const int C_MAX_SINGLE_LENGTH = 30;

const int D_MAX_HAIRPIN_LENGTH = 30;
const int D_MAX_BP_DIST_THRESHOLDS = 10;
const int D_MAX_BULGE_LENGTH = 30;
const int D_MAX_INTERNAL_LENGTH = 30;
const int D_MAX_INTERNAL_SYMMETRIC_LENGTH = 15;
const int D_MAX_INTERNAL_ASYMMETRY = 28;
const int D_MAX_INTERNAL_EXPLICIT_LENGTH = 4;
const int D_MAX_HELIX_LENGTH = 30;

const int BP_DIST_LAST_THRESHOLD = 132;
const int BP_DIST_THRESHOLDS[D_MAX_BP_DIST_THRESHOLDS] = { 3, 9, 12, 16, 21, 26, 34, 47, 71, BP_DIST_LAST_THRESHOLD };

#endif

// Local Variables:
// mode: C++
// c-basic-offset: 4
// End:
