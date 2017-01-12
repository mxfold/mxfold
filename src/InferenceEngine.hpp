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
private:
    const bool allow_noncomplementary;
    const int C_MAX_SINGLE_LENGTH;
    const int C_MIN_HAIRPIN_LENGTH;
    bool cache_initialized;
    const ParameterHash<RealT>* parameter_manager;
    ParameterHash<RealT>* parameter_count;

    //std::array<std::array<int, 256>, 256> is_complementary;
    
    // dimensions
    int L, SIZE;

    // sequence data
    std::vector<char> s;
    std::vector<int> offset;
    std::vector<int> allow_unpaired_position;
    std::vector<int> allow_unpaired, allow_paired;
    std::vector<RealT> loss_unpaired_position;
    std::vector<RealT> loss_unpaired, loss_paired;
    RealT loss_const;
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
    std::pair<RealT,RealT> cache_score_base_pair_dist[BP_DIST_LAST_THRESHOLD+1];
#endif
#if PARAMS_HAIRPIN_LENGTH
    std::pair<RealT,RealT> cache_score_hairpin_length[D_MAX_HAIRPIN_LENGTH+1];
#endif
#if PARAMS_HELIX_LENGTH
    std::pair<RealT,RealT> cache_score_helix_length[D_MAX_HELIX_LENGTH+1];
#endif

    // cache
    std::vector<std::vector<std::pair<RealT,RealT>>> cache_score_single;
    std::vector<std::pair<RealT,RealT> > cache_score_helix_sums;

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
    
    void CountUnpairedPosition(int i, RealT v);
    void CountUnpaired(int i,int j, RealT v);
    void CountIsolated(RealT v);
    void CountMultiBase(RealT v);
    void CountMultiPaired(RealT v);
    void CountMultiUnpaired(int i, RealT v);
    void CountExternalPaired(RealT v);
    void CountExternalUnpaired(int i, RealT v);
    void CountHelixStacking(int i,int j, RealT v);

    void CountJunctionA(int i, int j, RealT value);
    void CountJunctionB(int i, int j, RealT value);
    void CountBasePair(int i, int j, RealT value);
    void CountHairpin(int i, int j, RealT value);
    void CountHelix(int i, int j, int m, RealT value);
    void CountSingleNucleotides(int i, int j, int p, int q, RealT value);
    void CountSingle(int i, int j, int p, int q, RealT value);

    int EncodeTraceback(int i, int j) const;
    std::pair<int,int> DecodeTraceback(int s) const;

    void ClearCounts();
    void InitializeCache();
    void FinalizeCounts();
    
public:

    // constructor and destructor
    InferenceEngine(bool allow_noncomplementary, 
                    int max_single_length = DEFAULT_C_MAX_SINGLE_LENGTH,
                    int min_hairpin_length = DEFAULT_C_MIN_HAIRPIN_LENGTH);
    ~InferenceEngine();

    // load sequence
    void LoadSequence(const SStruct &sstruct);
    
    // load parameter values                        
    void LoadValues(const ParameterHash<RealT>* pm);
    
    // load loss function
    void UseLoss(const std::vector<int> &true_mapping, RealT example_loss);
    void UseLossBasePair(const std::vector<int> &true_mapping, RealT pos_w, RealT neg_w);
    void UseLossPosition(const std::vector<int> &true_mapping, RealT pos_w, RealT neg_w);
    void UseLossReactivity(const std::vector<float> &reactivity_unpair, const std::vector<float> &reactivity_pair, RealT pos_w, RealT neg_w);

    // use constraints
    void UseConstraints(const std::vector<int> &true_mapping);
    void UseSoftConstraints(const std::vector<float> &reactivity_unpair, const std::vector<float> &reactivity_pair, 
                            RealT threshold_unpaired_reactivity=0.7, RealT threshold_paired_reactivity=0.7, RealT scale_reactivity=0.1);

    // Viterbi inference
    void ComputeViterbi();
    RealT GetViterbiScore() const;
    std::vector<int> PredictPairingsViterbi() const;
    ParameterHash<RealT> ComputeViterbiFeatureCounts();

    // MEA inference
    void ComputeInside();
    RealT ComputeLogPartitionCoefficient() const;
    void ComputeOutside();
    ParameterHash<RealT> ComputeFeatureCountExpectations();
    void ComputePosterior();
    template <int GCE> std::vector<int> PredictPairingsPosterior(const RealT gamma) const;
    RealT *GetPosterior(const RealT posterior_cutoff) const;
};

#endif

// Local Variables:
// mode: C++
// c-basic-offset: 4
// End:
