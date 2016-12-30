//////////////////////////////////////////////////////////////////////
// InferenceEngine.hpp
//////////////////////////////////////////////////////////////////////

#ifndef INFERENCEENGINE_HPP
#define INFERENCEENGINE_HPP

#include <queue>
#include <vector>
#include <string>
#include <memory>
#include "Config.hpp"
#include "SStruct.hpp"
#include "ParameterHash.hpp"
#include "Utilities.hpp"
#include "LogSpace.hpp"
#include <iostream>
//////////////////////////////////////////////////////////////////////
// class InferenceEngine
//////////////////////////////////////////////////////////////////////

template<class RealT>
class InferenceEngine
{
public:
    typedef std::unique_ptr<ParameterHash<RealT>> ParamPtr;
    typedef std::unique_ptr<ParameterHash<uint>> CntPtr;

private:
    const bool allow_noncomplementary;
    std::array<std::array<char, 256>, 256> is_complementary;
    bool cache_initialized;
    ParamPtr parameter_manager;
    CntPtr parameter_count;
    
    // dimensions
    int L, SIZE;
#if PROFILE
    int N, SIZE2;
#endif

    // sequence data
    std::vector<char> s;
    std::vector<int> offset;
#if PROFILE
    std::vector<int> A;
    std::vector<RealT> weights;
#endif
    std::vector<int> allow_unpaired_position;
    std::vector<int> allow_unpaired, allow_paired;
    std::vector<RealT> loss_unpaired_position;
    std::vector<RealT> loss_unpaired, loss_paired;
    std::vector<float> reactivity_unpaired_position;
    std::vector<float> reactivity_unpaired, reactivity_paired;


    enum TRACEBACK_TYPE {
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
        TB_FN_HAIRPIN,
        TB_FN_SINGLE,
        TB_FN_BIFURCATION,
        TB_FE_STACKING,
        TB_FE_FN,
        TB_FC_FN,
        TB_FC_HELIX,
        TB_FC_FE,
#else
        TB_FC_HAIRPIN,
        TB_FC_SINGLE,
        TB_FC_BIFURCATION,
#endif
        TB_FM1_PAIRED,
        TB_FM1_UNPAIRED,
        TB_FM_BIFURCATION,
        TB_FM_UNPAIRED,
        TB_FM_FM1,
        TB_F5_ZERO,
        TB_F5_UNPAIRED,
        TB_F5_BIFURCATION,
        NUM_TRACEBACK_TYPES
    };
    
    // dynamic programming matrices
    std::vector<int> FCt, F5t, FMt, FM1t;            // traceback
    std::vector<RealT> FCv, F5v, FMv, FM1v;          // Viterbi  
    std::vector<RealT> FCi, F5i, FMi, FM1i;          // inside
    std::vector<RealT> FCo, F5o, FMo, FM1o;          // outside
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    std::vector<int> FEt, FNt;
    std::vector<RealT> FEv, FNv;
    std::vector<RealT> FEi, FNi;
    std::vector<RealT> FEo, FNo;
#endif
    
    std::vector<RealT> posterior;

    // cache
#if PARAMS_BASE_PAIR_DIST
    std::pair<RealT,uint> cache_score_base_pair_dist[BP_DIST_LAST_THRESHOLD+1];
#endif
#if PARAMS_HAIRPIN_LENGTH
    std::pair<RealT,uint> cache_score_hairpin_length[D_MAX_HAIRPIN_LENGTH+1];
#endif
#if PARAMS_HELIX_LENGTH
    std::pair<RealT,uint> cache_score_helix_length[D_MAX_HELIX_LENGTH+1];
#endif

#if PROFILE

    // multiple sequence scoring
#if PARAMS_BASE_PAIR
    std::vector<std::pair<RealT,RealT> > profile_score_base_pair;
#endif
#if PARAMS_TERMINAL_MISMATCH
    std::vector<std::pair<RealT,RealT> > profile_score_terminal_mismatch;
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_hairpin_3_nucleotides;
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_hairpin_4_nucleotides;
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_0x1_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_1x0_nucleotides;
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_0x2_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_2x0_nucleotides;
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_0x3_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_3x0_nucleotides;
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_1x1_nucleotides;
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_1x2_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_internal_2x1_nucleotides;
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_2x2_nucleotides;
#endif
#if PARAMS_HELIX_STACKING
    std::vector<std::pair<RealT,RealT> > profile_score_helix_stacking;
#endif
#if PARAMS_HELIX_CLOSING
    std::vector<std::pair<RealT,RealT> > profile_score_helix_closing;
#endif
#if PARAMS_DANGLE
    std::vector<std::pair<RealT,RealT> > profile_score_dangle_left;
    std::vector<std::pair<RealT,RealT> > profile_score_dangle_right;
#endif

    //自作パラメータ,multi
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_hairpin_5_nucleotides;
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_hairpin_6_nucleotides;
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_hairpin_7_nucleotides;
#endif

#if PARAMS_BULGE_0x4_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_0x4_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_4x0_nucleotides;
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_0x5_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_5x0_nucleotides;
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_0x6_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_bulge_6x0_nucleotides;
#endif

#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_1x3_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_internal_3x1_nucleotides;
#endif
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_2x3_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_internal_3x2_nucleotides;
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_3x3_nucleotides;
#endif

#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_1x4_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_internal_4x1_nucleotides;
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_2x4_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_internal_4x2_nucleotides;
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_3x4_nucleotides;
    std::vector<std::pair<RealT,RealT> > profile_score_internal_4x3_nucleotides;
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
    std::vector<std::pair<RealT,RealT> > profile_score_internal_4x4_nucleotides;
#endif
  
#endif

    // cache
    std::pair<RealT,uint> cache_score_single[C_MAX_SINGLE_LENGTH+1][C_MAX_SINGLE_LENGTH+1];
    std::vector<std::pair<RealT,uint> > cache_score_helix_sums;

    int ComputeRowOffset(int i, int N) const;
    bool IsComplementary(int i, int j) const;

    RealT ScoreUnpairedPosition(int i) const;
    RealT ScoreUnpaired(int i, int j) const;
    RealT ScoreIsolated() const;
    RealT ScoreMultiBase() const;
    RealT ScoreMultiPaired() const;
    RealT ScoreMultiUnpaired(int i) const;
    RealT ScoreExternalPaired() const;
    RealT ScoreExternalUnpaired(int i) const;
    RealT ScoreHelixStacking(int i, int j) const;

    RealT ScoreJunctionA(int i, int j) const;
    RealT ScoreJunctionB(int i, int j) const;
    RealT ScoreBasePair(int i, int j) const;
    RealT ScoreHairpin(int i, int j) const;
    RealT ScoreHelix(int i, int j, int m) const;
    RealT ScoreSingleNucleotides(int i, int j, int p, int q) const;
    RealT ScoreSingle(int i, int j, int p, int q) const;
    
    void CountUnpairedPosition(int i, uint v);
    void CountUnpaired(int i,int j, uint v);
    void CountIsolated(uint v);
    void CountMultiBase(uint v);
    void CountMultiPaired(uint v);
    void CountMultiUnpaired(int i, uint v);
    void CountExternalPaired(uint v);
    void CountExternalUnpaired(int i, uint v);
    void CountHelixStacking(int i,int j, uint v);

    void CountJunctionA(int i, int j, uint value);
    void CountJunctionB(int i, int j, uint value);
    void CountBasePair(int i, int j, uint value);
    void CountHairpin(int i, int j, uint value);
    void CountHelix(int i, int j, int m, uint value);
    void CountSingleNucleotides(int i, int j, int p, int q, uint value);
    void CountSingle(int i, int j, int p, int q, uint value);

    int EncodeTraceback(int i, int j) const;
    std::pair<int,int> DecodeTraceback(int s) const;

    void ClearCounts();
    void InitializeCache();
    void FinalizeCounts();
    
#if PROFILE
    void ComputeProfileScore(RealT &profile_score, const int *pos, int dimensions, std::pair<RealT,RealT> *table);
    void ConvertProfileCount(const RealT &profile_score, const int *pos, int dimensions, std::pair<RealT,RealT> *table);
#endif
    
public:

    // constructor and destructor
    InferenceEngine(bool allow_noncomplementary);
    ~InferenceEngine();

    // load sequence
    void LoadSequence(const SStruct &sstruct, int use_reactivity=0);
    
    // load parameter values                        
    ParamPtr LoadValues(ParamPtr pm);
    
    // load loss function
    void UseLoss(const std::vector<int> &true_mapping, RealT example_loss);

    // use constraints
    void UseConstraints(const std::vector<int> &true_mapping);

    // Viterbi inference
    void ComputeViterbi();
    RealT GetViterbiScore() const;
    std::vector<int> PredictPairingsViterbi() const;
    CntPtr ComputeViterbiFeatureCounts();

    // MEA inference
    void ComputeInside();
    RealT ComputeLogPartitionCoefficient() const;
    void ComputeOutside();
    CntPtr ComputeFeatureCountExpectations();
    void ComputePosterior();
    std::vector<int> PredictPairingsPosterior(const RealT gamma) const;
    RealT *GetPosterior(const RealT posterior_cutoff) const;
};

#endif

// Local Variables:
// mode: C++
// c-basic-offset: 4
// End:
