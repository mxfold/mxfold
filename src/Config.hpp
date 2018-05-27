//////////////////////////////////////////////////////////////////////
// Config.hpp
//
// Global configuration file.
//////////////////////////////////////////////////////////////////////

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <vector>

typedef float param_value_type;

typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<VVI> VVVI;
typedef std::vector<VVVI> VVVVI;

//////////////////////////////////////////////////////////////////////
// Options related to general inference
//////////////////////////////////////////////////////////////////////

// showing timings for inference routines
#define SHOW_TIMINGS                               0

// use straightforward calculation for FM2 matrix
#define SIMPLE_FM2                                 0

// use candidate list optimization for Viterbi parsing
#define CANDIDATE_LIST                             1

// use caching algorithm for fast helix length scores
#define FAST_HELIX_LENGTHS                         1

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

#else

// our full model
#define PARAMS_BASE_PAIR                           1
#define PARAMS_BASE_PAIR_DIST                      1
#define PARAMS_TERMINAL_MISMATCH                   1
#define PARAMS_HAIRPIN_LENGTH                      1
#define PARAMS_HELIX_LENGTH                        1
#define PARAMS_ISOLATED_BASE_PAIR                  1
#define PARAMS_INTERNAL_EXPLICIT                   1
#define PARAMS_BULGE_LENGTH                        1
#define PARAMS_INTERNAL_LENGTH                     1
#define PARAMS_INTERNAL_SYMMETRY                   1
#define PARAMS_INTERNAL_ASYMMETRY                  1
#define PARAMS_HELIX_STACKING                      1
#define PARAMS_HELIX_CLOSING                       1
#define PARAMS_MULTI_LENGTH                        1
#define PARAMS_DANGLE                              1
#define PARAMS_EXTERNAL_LENGTH                     1

// extended parameters
#define PARAMS_INTERNAL_NUCLEOTIDES                1
#define PARAMS_HAIRPIN_NUCLEOTIDES                 1
#define PARAMS_TERMINAL_MISMATCH_HAIRPIN           0
#define PARAMS_TERMINAL_MISMATCH_INTERNAL          0
#define PARAMS_TERMINAL_MISMATCH_INTERNAL_1N       0
#define PARAMS_TERMINAL_MISMATCH_INTERNAL_23       0
#define PARAMS_TERMINAL_MISMATCH_MULTI             0
#define PARAMS_TERMINAL_MISMATCH_EXTERNAL          0

//#define PARAMS_VIENNA_COMPAT                     1

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

const int DEFAULT_C_MIN_HAIRPIN_LENGTH = 0;
const int DEFAULT_C_MIN_HAIRPIN_LENGTH_PREDICT = 3;
const int DEFAULT_C_MAX_HAIRPIN_NUCLEOTIDES_LENGTH = 7;
const int DEFAULT_C_MAX_SINGLE_LENGTH = 40;
const int DEFAULT_C_MAX_SINGLE_NUCLEOTIDES_LENGTH = 7;

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
