///////////////////////////////////////////////////////////////////////
// InferenceEngine.cpp
//////////////////////////////////////////////////////////////////////

#include "InferenceEngine.hpp"

//////////////////////////////////////////////////////////////////////
// UPDATE_MAX()
//
// Macro for updating a score/traceback pointer which does not
// evaluate t unless an update is needed.  Make sure that this is
// used as a stand-alone statement (i.e., not the "if" condition
// of an if-then-else statement.)
//////////////////////////////////////////////////////////////////////

#define UPDATE_MAX(bs,bt,s,t) { RealT work(s); if ((work)>(bs)) { (bs)=(work); (bt)=(t); } }

//////////////////////////////////////////////////////////////////////
// FillScores()
// FillCounts()
// 
// Routines for setting scores and counts quickly. 
//////////////////////////////////////////////////////////////////////

template<class ITR, class V>
void FillScores(ITR begin, ITR end, V value)
{
    while (begin != end)
    {
        begin->first = value;
        ++begin;
    }
}

template<class ITR, class V>
void FillCounts(ITR begin, ITR end, V value)
{
    while (begin != end)
    {
        begin->second = value;
        ++begin;
    }
}

//////////////////////////////////////////////////////////////////////
// ComputeRowOffset()
//
// Consider an N x N upper triangular matrix whose elements are
// stored in a one-dimensional flat array using the following
// row-major indexing scheme:
//
//     0  1  2  3     <-- row 0
//        4  5  6     <-- row 1
//           7 [8]    <-- row 2
//              9     <-- row 3
//
// Assuming 0-based indexing, this function computes offset[i]
// for the ith row such that offset[i]+j is the index of the
// (i,j)th element of the upper triangular matrix in the flat
// array.
//
// For example, offset[2] = 5, so the (2,3)th element of the
// upper triangular matrix (marked in the picture above) can be 
// found at position offset[2]+3 = 5+3 = 8 in the flat array.
//////////////////////////////////////////////////////////////////////

template<class RealT>
int InferenceEngine<RealT>::ComputeRowOffset(int i, int N) const
{
    Assert(i >= 0 && i <= N, "Index out-of-bounds.");
    
    // equivalent to:
    //   return N*(N+1)/2 - (N-i)*(N-i+1)/2 - i;
    return i*(N+N-i-1)/2;
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::IsComplementary()
//
// Determine if a pair of positions is considered "complementary."
//////////////////////////////////////////////////////////////////////

template<class RealT>
bool InferenceEngine<RealT>::IsComplementary(int i, int j) const
{
    Assert(1 <= i && i <= L, "Index out-of-bounds.");
    Assert(1 <= j && j <= L, "Index out-of-bounds.");

#if !PROFILE
    return is_complementary[s[i]][s[j]];
#else
    RealT complementary_weight = 0;
    RealT total_weight = 0;

    for (int k = 0; k < N; k++)
    {
        if (is_complementary[A[k*(L+1)+i]][A[k*(L+1)+j]]) complementary_weight += weights[k];
        total_weight += weights[k];
    }

    return complementary_weight / total_weight >= std::min(RealT(N-1) / RealT(N), RealT(0.5));
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::InferenceEngine()
//
// Constructor
//////////////////////////////////////////////////////////////////////

template<class RealT>
InferenceEngine<RealT>::InferenceEngine(bool allow_noncomplementary) :
    allow_noncomplementary(allow_noncomplementary),
    cache_initialized(false),
    L(0),
    SIZE(0)
#if PROFILE
    , N(0)
    , SIZE2(0)
#endif

{
    // precompute complementary pairings
    for (auto e : is_complementary)
        std::fill(std::begin(e), std::end(e), 0);
    
    is_complementary['A']['U'] = is_complementary['U']['A'] = 1;
    is_complementary['G']['U'] = is_complementary['U']['G'] = 1;
    is_complementary['C']['G'] = is_complementary['G']['C'] = 1;
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::~InferenceEngine()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
InferenceEngine<RealT>::~InferenceEngine()
{}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::LoadSequence()
//
// Load an RNA sequence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::LoadSequence(const SStruct &sstruct, int use_reactivity)
{
    //std::cout <<"reactivity in InferenceEngine.ipp is " << use_reactivity << std::endl; 
    const std::vector<float> &reactivity_unpair= sstruct.GetReactivityUnpair();
    cache_initialized = false;
    
    // compute dimensions
    L = sstruct.GetLength();
    SIZE = (L+1)*(L+2) / 2;
#if PROFILE
    N = sstruct.GetNumSequences();
    SIZE2 = (L+1)*(L+1);
#endif
    
    // allocate memory
    s.resize(L+1);
#if PROFILE
    A.resize(N*(L+1));
    weights.resize(N);
#endif
    offset.resize(L+1);
    allow_unpaired_position.resize(L+1);
    allow_unpaired.resize(SIZE);
    allow_paired.resize(SIZE);
    loss_unpaired_position.resize(L+1);
    loss_unpaired.resize(SIZE);
    loss_paired.resize(SIZE);
    
    reactivity_unpaired_position.resize(L+1);
    reactivity_unpaired.resize(SIZE);
    reactivity_paired.resize(SIZE);
    
#if PROFILE

#if PARAMS_BASE_PAIR
    profile_score_base_pair.clear();                 profile_score_base_pair.resize(SIZE2);
#endif
#if PARAMS_TERMINAL_MISMATCH
    profile_score_terminal_mismatch.clear();         profile_score_terminal_mismatch.resize(SIZE2);
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
    profile_score_hairpin_3_nucleotides.clear();     profile_score_hairpin_3_nucleotides.resize(L+1);
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
    profile_score_hairpin_4_nucleotides.clear();     profile_score_hairpin_4_nucleotides.resize(L+1);
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
    profile_score_bulge_0x1_nucleotides.clear();     profile_score_bulge_0x1_nucleotides.resize(L+1);
    profile_score_bulge_1x0_nucleotides.clear();     profile_score_bulge_1x0_nucleotides.resize(L+1);
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
    profile_score_bulge_0x2_nucleotides.clear();     profile_score_bulge_0x2_nucleotides.resize(L+1);
    profile_score_bulge_2x0_nucleotides.clear();     profile_score_bulge_2x0_nucleotides.resize(L+1);
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
    profile_score_bulge_0x3_nucleotides.clear();     profile_score_bulge_0x3_nucleotides.resize(L+1);
    profile_score_bulge_3x0_nucleotides.clear();     profile_score_bulge_3x0_nucleotides.resize(L+1);
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
    profile_score_internal_1x1_nucleotides.clear();  profile_score_internal_1x1_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
    profile_score_internal_1x2_nucleotides.clear();  profile_score_internal_1x2_nucleotides.resize(SIZE2);
    profile_score_internal_2x1_nucleotides.clear();  profile_score_internal_2x1_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
    profile_score_internal_2x2_nucleotides.clear();  profile_score_internal_2x2_nucleotides.resize(SIZE2);
#endif
#if PARAMS_HELIX_STACKING
    profile_score_helix_stacking.clear();            profile_score_helix_stacking.resize(SIZE2);
#endif
#if PARAMS_HELIX_CLOSING
    profile_score_helix_closing.clear();             profile_score_helix_closing.resize(SIZE2);
#endif
#if PARAMS_DANGLE
    profile_score_dangle_left.clear();               profile_score_dangle_left.resize(SIZE2);
    profile_score_dangle_right.clear();              profile_score_dangle_right.resize(SIZE2);
#endif
    //自作パラメータ
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
    profile_score_hairpin_5_nucleotides.clear();     profile_score_hairpin_5_nucleotides.resize(L+1);
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
    profile_score_hairpin_6_nucleotides.clear();     profile_score_hairpin_6_nucleotides.resize(L+1);
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
    profile_score_hairpin_7_nucleotides.clear();     profile_score_hairpin_7_nucleotides.resize(L+1);
#endif

#if PARAMS_BULGE_0x4_NUCLEOTIDES
    profile_score_bulge_0x4_nucleotides.clear();     profile_score_bulge_0x4_nucleotides.resize(L+1);
    profile_score_bulge_4x0_nucleotides.clear();     profile_score_bulge_4x0_nucleotides.resize(L+1);
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
    profile_score_bulge_0x5_nucleotides.clear();     profile_score_bulge_0x5_nucleotides.resize(L+1);
    profile_score_bulge_5x0_nucleotides.clear();     profile_score_bulge_5x0_nucleotides.resize(L+1);
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
    profile_score_bulge_0x6_nucleotides.clear();     profile_score_bulge_0x6_nucleotides.resize(L+1);
    profile_score_bulge_6x0_nucleotides.clear();     profile_score_bulge_6x0_nucleotides.resize(L+1);
#endif

#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
    profile_score_internal_1x3_nucleotides.clear();  profile_score_internal_1x3_nucleotides.resize(SIZE2);
    profile_score_internal_3x1_nucleotides.clear();  profile_score_internal_3x1_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
    profile_score_internal_2x3_nucleotides.clear();  profile_score_internal_2x3_nucleotides.resize(SIZE2);
    profile_score_internal_3x2_nucleotides.clear();  profile_score_internal_3x2_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
    profile_score_internal_3x3_nucleotides.clear();  profile_score_internal_3x3_nucleotides.resize(SIZE2);
#endif
 
#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
    profile_score_internal_1x4_nucleotides.clear();  profile_score_internal_1x4_nucleotides.resize(SIZE2);
    profile_score_internal_4x1_nucleotides.clear();  profile_score_internal_4x1_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
    profile_score_internal_2x4_nucleotides.clear();  profile_score_internal_2x4_nucleotides.resize(SIZE2);
    profile_score_internal_4x2_nucleotides.clear();  profile_score_internal_4x2_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
    profile_score_internal_3x4_nucleotides.clear();  profile_score_internal_3x4_nucleotides.resize(SIZE2);
    profile_score_internal_4x3_nucleotides.clear();  profile_score_internal_4x3_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
    profile_score_internal_4x4_nucleotides.clear();  profile_score_internal_4x4_nucleotides.resize(SIZE2);
#endif
    
#endif

#if FAST_HELIX_LENGTHS
    cache_score_helix_sums.clear();                  cache_score_helix_sums.resize((2*L+1)*L);
#endif

    // convert sequences to index representation
    const std::string &sequence = sstruct.GetSequences()[0];
    s[0] = -1; //BYTE(alphabet.size());
    for (int i = 1; i <= L; i++)
    {
        s[i] = toupper(sequence[i]);
    }
#if PROFILE
    const std::vector<std::string> &alignment = sstruct.GetSequences();
    for (int k = 0; k < N; k++)
    {
        A[k*(L+1)+0] = -1;
        for (int i = 1; i <= L; i++)
        {
            A[k*(L+1)+i] = alignment[k][i];
        }
    }

    weights = ConvertVector<RealT>(sstruct.ComputePositionBasedSequenceWeights());
#endif
    
    // compute indexing scheme for upper triangular arrays;
    // also allow each position to be unpaired by default, and
    // set the loss for each unpaired position to zero
    for (int i = 0; i <= L; i++)
    {
        offset[i] = ComputeRowOffset(i,L+1);
        allow_unpaired_position[i] = 1;
        loss_unpaired_position[i] = RealT(0);
        //std::cout << reactivity_unpair[i] << std::endl;
        //std::cout <<"reactivity in InferenceEngine.ipp is " << use_reactivity << std::endl; 
        if(reactivity_unpair[i]>0.7){
            reactivity_unpaired_position[i] =  0.1*use_reactivity * reactivity_unpair[i];
            //std::cout << reactivity_unpaired_position[i] << std::endl;
        }
        else
        {
            reactivity_unpaired_position[i] =  0;
        }
    }

    //std::cout << reactivity_unpaired_position[18] << std::endl;
    for (int i = 0; i <= L; i++)
    {
        reactivity_unpaired[offset[i]+i] = RealT(0);
        reactivity_paired[offset[i]+i] = RealT(NEG_INF);
        for (int j = i+1; j <= L; j++)
        {
            reactivity_unpaired[offset[i]+j] = 
                reactivity_unpaired[offset[i]+j-1] +
                reactivity_unpaired_position[j];
            reactivity_paired[offset[i]+j] =0; 
        }
    }

    // allow all ranges to be unpaired, and all pairs of letters
    // to be paired; set the respective losses to zero    
    for (int i = 0; i < SIZE; i++)
    {
        allow_unpaired[i] = 1;
        allow_paired[i] = 1;
        loss_unpaired[i] = RealT(0);
        loss_paired[i] = RealT(0);
    }

    // prevent the non-letter before each sequence from pairing with anything;
    // also prevent each letter from pairing with itself
    for (int i = 0; i <= L; i++)
    {
        allow_paired[offset[0]+i] = 0;
        allow_paired[offset[i]+i] = 0;
    }

    // enforce complementarity of base-pairings
    if (!allow_noncomplementary)
    {
        // for each pair of non-complementary letters in the sequence, disallow the pairing
        for (int i = 1; i <= L; i++)
        {
            for (int j = i+1; j <= L; j++)
            {
                if (!IsComplementary(i,j))
                    allow_paired[offset[i]+j] = 0;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::InitializeCache()
//
// Initialize scoring cache prior to inference.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::InitializeCache()
{
    if (cache_initialized) return;
    cache_initialized = true;
    const auto& pm = *parameter_manager;

    // initialize length and distance scoring
#if PARAMS_BASE_PAIR_DIST
    for (int j = 0; j <= BP_DIST_LAST_THRESHOLD; j++)
        cache_score_base_pair_dist[j].first = RealT(0);
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
        for (int j = BP_DIST_THRESHOLDS[i]; j <= BP_DIST_LAST_THRESHOLD; j++)
            cache_score_base_pair_dist[j].first += pm.base_pair_dist_at_least(i);
#endif
    
#if PARAMS_HAIRPIN_LENGTH
    cache_score_hairpin_length[0].first = pm.hairpin_length_at_least(0);
    for (int i = 1; i <= D_MAX_HAIRPIN_LENGTH; i++)
        cache_score_hairpin_length[i].first = cache_score_hairpin_length[i-1].first + pm.hairpin_length_at_least(i);
#endif

#if PARAMS_HELIX_LENGTH
    cache_score_helix_length[0].first = pm.helix_length_at_least(0);
    for (int i = 1; i <= D_MAX_HELIX_LENGTH; i++)
        cache_score_helix_length[i].first = cache_score_helix_length[i-1].first + pm.helix_length_at_least(i);
#endif

#if PARAMS_BULGE_LENGTH
    RealT temp_cache_score_bulge_length[D_MAX_BULGE_LENGTH+1];
    temp_cache_score_bulge_length[0] = pm.bulge_length_at_least(0);
    for (int i = 1; i <= D_MAX_BULGE_LENGTH; i++)
        temp_cache_score_bulge_length[i] = temp_cache_score_bulge_length[i-1] + pm.bulge_length_at_least(i);
#endif
    
#if PARAMS_INTERNAL_LENGTH
    RealT temp_cache_score_internal_length[D_MAX_INTERNAL_LENGTH+1];
    temp_cache_score_internal_length[0] = pm.internal_length_at_least(0);
    for (int i = 1; i <= D_MAX_INTERNAL_LENGTH; i++)
        temp_cache_score_internal_length[i] = temp_cache_score_internal_length[i-1] + pm.internal_length_at_least(i);
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    RealT temp_cache_score_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
    temp_cache_score_internal_symmetric_length[0] = pm.internal_symmetric_length_at_least(0);
    for (int i = 1; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
        temp_cache_score_internal_symmetric_length[i] = temp_cache_score_internal_symmetric_length[i-1] + pm.internal_symmetric_length_at_least(i);
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    RealT temp_cache_score_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
    temp_cache_score_internal_asymmetry[0] = pm.internal_asymmetry_at_least(0);
    for (int i = 1; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
        temp_cache_score_internal_asymmetry[i] = temp_cache_score_internal_asymmetry[i-1] + pm.internal_asymmetry_at_least(i);
#endif
    
    // precompute score for single-branch loops of length l1 and l2
    //長さl1とl2の内部ループスコア
    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
    {
        for (int l2 = 0; l1+l2 <= C_MAX_SINGLE_LENGTH; l2++)
        {
            cache_score_single[l1][l2].first = RealT(0);

            // skip over stacking pairs
            if (l1 == 0 && l2 == 0) continue;

            // consider bulge loops
            if (l1 == 0 || l2 == 0)
            {
#if PARAMS_BULGE_LENGTH
                cache_score_single[l1][l2].first += temp_cache_score_bulge_length[std::min(D_MAX_BULGE_LENGTH, l1+l2)];
#endif
            }

            // consider internal loops
            else
            {
#if PARAMS_INTERNAL_EXPLICIT
                if (l1 <= D_MAX_INTERNAL_EXPLICIT_LENGTH && l2 <= D_MAX_INTERNAL_EXPLICIT_LENGTH)
                    cache_score_single[l1][l2].first += pm.internal_explicit(l1, l2);
#endif
#if PARAMS_INTERNAL_LENGTH
                cache_score_single[l1][l2].first += temp_cache_score_internal_length[std::min(D_MAX_INTERNAL_LENGTH, l1+l2)];
#endif
#if PARAMS_INTERNAL_SYMMETRY
                if (l1 == l2)
                    cache_score_single[l1][l2].first += temp_cache_score_internal_symmetric_length[std::min(D_MAX_INTERNAL_SYMMETRIC_LENGTH, l1)];
#endif
#if PARAMS_INTERNAL_ASYMMETRY
                cache_score_single[l1][l2].first += temp_cache_score_internal_asymmetry[std::min(D_MAX_INTERNAL_ASYMMETRY, Abs(l1-l2))];
#endif
            }
        }
    }
    
#if PROFILE
    // initialize counts for profile scoring
    for (int i = 0; i <= L; i++)
    {
        for (int j = 0; j <= L; j++)
        {
#if PARAMS_BASE_PAIR
            {
                const int pos[2] = {i, j};
                ComputeProfileScore(profile_score_base_pair[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_base_pair));
            }
#endif
#if PARAMS_TERMINAL_MISMATCH
            {
                const int pos[4] = {i, j+1, i+1, j};
                ComputeProfileScore(profile_score_terminal_mismatch[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_terminal_mismatch));
            }
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ComputeProfileScore(profile_score_hairpin_3_nucleotides[i].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_3_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[4] = {i+1, i+2, i+3, i+4};
                ComputeProfileScore(profile_score_hairpin_4_nucleotides[i].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_4_nucleotides));
            }
#endif
	    //自作パラメータ
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[5] = {i+1, i+2, i+3, i+4,i+5};
                ComputeProfileScore(profile_score_hairpin_5_nucleotides[i].first, pos, 5, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_5_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[6] = {i+1, i+2, i+3, i+4,i+5,i+6};
                ComputeProfileScore(profile_score_hairpin_6_nucleotides[i].first, pos, 6, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_6_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[7] = {i+1, i+2, i+3, i+4,i+5,i+6,i+7};
                ComputeProfileScore(profile_score_hairpin_7_nucleotides[i].first, pos, 7, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_7_nucleotides));
            }
#endif
	    //自作終わり

#if PARAMS_BULGE_0x1_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[1] = {j};
                ComputeProfileScore(profile_score_bulge_0x1_nucleotides[j].first, pos, 1, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x1_nucleotides));
            }
            if (j == 0)
            {
                const int pos[1] = {i+1};
                ComputeProfileScore(profile_score_bulge_1x0_nucleotides[i].first, pos, 1, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_1x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[2] = {j-1, j};
                ComputeProfileScore(profile_score_bulge_0x2_nucleotides[j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x2_nucleotides));
            }
            if (j == 0)
            {
                const int pos[2] = {i+1, i+2};
                ComputeProfileScore(profile_score_bulge_2x0_nucleotides[i].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_2x0_nucleotides));
            }
#endif            
#if PARAMS_BULGE_0x3_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[3] = {j-2, j-1, j};
                ComputeProfileScore(profile_score_bulge_0x3_nucleotides[j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x3_nucleotides));
            }
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ComputeProfileScore(profile_score_bulge_3x0_nucleotides[i].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_3x0_nucleotides));
            }
#endif
	    //自作パラメータ
#if PARAMS_BULGE_0x4_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[4] = {j-3,j-2, j-1, j};
                ComputeProfileScore(profile_score_bulge_0x4_nucleotides[j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x4_nucleotides));
            }
            if (j == 0)
            {
                const int pos[4] = {i+1, i+2, i+3,i+4};
                ComputeProfileScore(profile_score_bulge_4x0_nucleotides[i].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_4x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[5] = {j-4,j-3,j-2, j-1, j};
                ComputeProfileScore(profile_score_bulge_0x5_nucleotides[j].first, pos, 5, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x5_nucleotides));
            }
            if (j == 0)
            {
                const int pos[5] = {i+1, i+2, i+3,i+4,i+5};
                ComputeProfileScore(profile_score_bulge_5x0_nucleotides[i].first, pos, 5, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_5x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[6] = {j-5,j-4,j-3,j-2, j-1, j};
                ComputeProfileScore(profile_score_bulge_0x6_nucleotides[j].first, pos, 6, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x6_nucleotides));
            }
            if (j == 0)
            {
                const int pos[6] = {i+1, i+2, i+3,i+4,i+5,i+6};
                ComputeProfileScore(profile_score_bulge_6x0_nucleotides[i].first, pos, 6, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_6x0_nucleotides));
            }
#endif
	    //自作終了            
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
            {
                const int pos[2] = {i+1, j};
                ComputeProfileScore(profile_score_internal_1x1_nucleotides[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_1x1_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
            {
                const int pos[3] = {i+1, j-1, j};
                ComputeProfileScore(profile_score_internal_1x2_nucleotides[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_1x2_nucleotides));
            }
            {
                const int pos[3] = {i+1, i+2, j};
                ComputeProfileScore(profile_score_internal_2x1_nucleotides[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_2x1_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
            {
                const int pos[4] = {i+1, i+2, j-1, j};
                ComputeProfileScore(profile_score_internal_2x2_nucleotides[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_2x2_nucleotides));
            }
#endif
	    //自作パラメータ
#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
            {
                const int pos[4] = {i+1, j-2, j-1, j};
                ComputeProfileScore(profile_score_internal_1x3_nucleotides[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_1x3_nucleotides));
            }
            {
                const int pos[4] = {i+1, i+2, i+3, j};
                ComputeProfileScore(profile_score_internal_3x1_nucleotides[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_3x1_nucleotides));
            }
#endif     
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
            {
                const int pos[5] = {i+1, i+2, j-2, j-1, j};
                ComputeProfileScore(profile_score_internal_2x3_nucleotides[i*(L+1)+j].first, pos, 5, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_2x3_nucleotides));
            }
            {
                const int pos[5] = {i+1, i+2, i+3, j-1, j};
                ComputeProfileScore(profile_score_internal_3x2_nucleotides[i*(L+1)+j].first, pos, 5, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_3x2_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
            {
                const int pos[6] = {i+1, i+2, i+3, j-2, j-1, j};
                ComputeProfileScore(profile_score_internal_3x3_nucleotides[i*(L+1)+j].first, pos, 6, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_3x3_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
	    {
                const int pos[5] = {i+1, j-3, j-2, j-1, j};
                ComputeProfileScore(profile_score_internal_1x4_nucleotides[i*(L+1)+j].first, pos, 5, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_1x4_nucleotides));
	    }
	    {
                const int pos[5] = {i+1, i+2, i+3, i+4, j};
                ComputeProfileScore(profile_score_internal_4x1_nucleotides[i*(L+1)+j].first, pos, 5, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_4x1_nucleotides));
	    }
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
	    {
                const int pos[6] = {i+1, i+2, j-3, j-2, j-1, j};
                ComputeProfileScore(profile_score_internal_2x4_nucleotides[i*(L+1)+j].first, pos, 6, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_2x4_nucleotides));
	    }
	    {
                const int pos[6] = {i+1, i+2, i+3, i+4, j-1, j};
                ComputeProfileScore(profile_score_internal_4x2_nucleotides[i*(L+1)+j].first, pos, 6, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_4x2_nucleotides));
	    }
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
	    {
                const int pos[7] = {i+1, i+2, i+3, j-3, j-2, j-1, j};
                ComputeProfileScore(profile_score_internal_3x4_nucleotides[i*(L+1)+j].first, pos, 7, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_3x4_nucleotides));
	    }
	    {
                const int pos[7] = {i+1, i+2, i+3, i+4, j-2, j-1, j};
                ComputeProfileScore(profile_score_internal_4x3_nucleotides[i*(L+1)+j].first, pos, 7, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_4x3_nucleotides));
	    }
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
	    {
                const int pos[8] = {i+1, i+2, i+3, i+4, j-3, j-2, j-1, j};
                ComputeProfileScore(profile_score_internal_4x4_nucleotides[i*(L+1)+j].first, pos, 8, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_4x4_nucleotides));
	    }
#endif
	    //自作終了

#if PARAMS_HELIX_STACKING
            {
                const int pos[4] = {i, j, i+1, j-1};
                ComputeProfileScore(profile_score_helix_stacking[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_helix_stacking));
            }
#endif
#if PARAMS_HELIX_CLOSING
            {
                const int pos[2] = {i, j+1};
                ComputeProfileScore(profile_score_helix_closing[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_helix_closing));
            }
#endif
#if PARAMS_DANGLE
            {
                const int pos[3] = {i, j+1, i+1};
                ComputeProfileScore(profile_score_dangle_left[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_dangle_left));
            }
            {
                const int pos[3] = {i, j+1, j};
                ComputeProfileScore(profile_score_dangle_right[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_dangle_right));
            }
#endif
        }
    }

#endif

#if FAST_HELIX_LENGTHS
    // precompute helix partial sums
    FillScores(cache_score_helix_sums.begin(), cache_score_helix_sums.end(), 0);
    for (int i = L; i >= 1; i--)
    {
        for (int j = i+3; j <= L; j++)
        {
            cache_score_helix_sums[(i+j)*L+j-i].first = cache_score_helix_sums[(i+j)*L+j-i-2].first;
            if (allow_paired[offset[i+1]+j-1])
            {
                cache_score_helix_sums[(i+j)*L+j-i].first += ScoreBasePair(i+1,j-1);
                if (allow_paired[offset[i]+j])
                    cache_score_helix_sums[(i+j)*L+j-i].first += ScoreHelixStacking(i,j);
            }
        }
    }
#endif

}

#if PROFILE

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeProfileScore()
//
// Compute profile score for a single location.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputeProfileScore(RealT &profile_score, const int *pos, int dimensions, std::pair<RealT,RealT> *table)
{
    profile_score = 0;

    // consider all sequences
    for (int k = 0; k < N; k++)
    {
        bool valid = true;
        int index = 0;
        int *seq = &A[k*(L+1)];
        
        // extract letters of the pattern for the current sequence
        for (int d = 0; valid && d < dimensions; d++)
        {
            if (pos[d] < 1 || pos[d] > L)
                valid = false;
            else
            {
                BYTE c = seq[pos[d]];
                if (c == BYTE(alphabet.size()))
                    valid = false;
                else
                    index = index * (M+1) + c;
            }
        }

        // add contribution of pattern to score
        if (valid) profile_score += weights[k] * table[index].first;
    }
}

#endif
            
//////////////////////////////////////////////////////////////////////
// InferenceEngine::LoadValues()
//
// Load parameter values.
//////////////////////////////////////////////////////////////////////

template<class RealT>
typename InferenceEngine<RealT>::ParamPtr
InferenceEngine<RealT>::LoadValues(ParamPtr pm)
{
    cache_initialized = false;
    std::swap(pm, parameter_manager);
    return std::move(pm);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ClearCounts()
//
// Set all counts to zero.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ClearCounts()
{
    // clear counts for cache
#if PARAMS_BASE_PAIR_DIST
    for (int i = 0; i <= BP_DIST_LAST_THRESHOLD; i++)
        cache_score_base_pair_dist[i].second = 0;
#endif

#if PARAMS_HAIRPIN_LENGTH
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
        cache_score_hairpin_length[i].second = 0;
#endif

#if PARAMS_HELIX_LENGTH
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
        cache_score_helix_length[i].second = 0;
#endif

    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
        for (int l2 = 0; l2 <= C_MAX_SINGLE_LENGTH; l2++)
            cache_score_single[l1][l2].second = 0;
    
#if FAST_HELIX_LENGTHS
    FillCounts(cache_score_helix_sums.begin(), cache_score_helix_sums.end(), 0);
#endif

    // clear counts for profiles
#if PROFILE

#if PARAMS_BASE_PAIR
    FillCounts(profile_score_base_pair.begin(), profile_score_base_pair.end(), RealT(0));
#endif
#if PARAMS_TERMINAL_MISMATCH
    FillCounts(profile_score_terminal_mismatch.begin(), profile_score_terminal_mismatch.end(), RealT(0));
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
    FillCounts(profile_score_hairpin_3_nucleotides.begin(), profile_score_hairpin_3_nucleotides.end(), RealT(0));
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
    FillCounts(profile_score_hairpin_4_nucleotides.begin(), profile_score_hairpin_4_nucleotides.end(), RealT(0));
#endif
    //自作パラメータ
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
    FillCounts(profile_score_hairpin_5_nucleotides.begin(), profile_score_hairpin_5_nucleotides.end(), RealT(0));
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
    FillCounts(profile_score_hairpin_6_nucleotides.begin(), profile_score_hairpin_6_nucleotides.end(), RealT(0));
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
    FillCounts(profile_score_hairpin_7_nucleotides.begin(), profile_score_hairpin_7_nucleotides.end(), RealT(0));
#endif
    //自作終わり
#if PARAMS_BULGE_0x1_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x1_nucleotides.begin(), profile_score_bulge_0x1_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_1x0_nucleotides.begin(), profile_score_bulge_1x0_nucleotides.end(), RealT(0));
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x2_nucleotides.begin(), profile_score_bulge_0x2_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_2x0_nucleotides.begin(), profile_score_bulge_2x0_nucleotides.end(), RealT(0));
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x3_nucleotides.begin(), profile_score_bulge_0x3_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_3x0_nucleotides.begin(), profile_score_bulge_3x0_nucleotides.end(), RealT(0));
#endif
    //自作パラメータ
#if PARAMS_BULGE_0x4_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x4_nucleotides.begin(), profile_score_bulge_0x4_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_4x0_nucleotides.begin(), profile_score_bulge_4x0_nucleotides.end(), RealT(0));
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x5_nucleotides.begin(), profile_score_bulge_0x5_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_5x0_nucleotides.begin(), profile_score_bulge_5x0_nucleotides.end(), RealT(0));
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x6_nucleotides.begin(), profile_score_bulge_0x6_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_6x0_nucleotides.begin(), profile_score_bulge_6x0_nucleotides.end(), RealT(0));
#endif
    //自作終わり
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
    FillCounts(profile_score_internal_1x1_nucleotides.begin(), profile_score_internal_1x1_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
    FillCounts(profile_score_internal_1x2_nucleotides.begin(), profile_score_internal_1x2_nucleotides.end(), RealT(0));
    FillCounts(profile_score_internal_2x1_nucleotides.begin(), profile_score_internal_2x1_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
    FillCounts(profile_score_internal_2x2_nucleotides.begin(), profile_score_internal_2x2_nucleotides.end(), RealT(0));
#endif
    //自作パラメータ
#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
    FillCounts(profile_score_internal_1x3_nucleotides.begin(), profile_score_internal_1x3_nucleotides.end(), RealT(0));
    FillCounts(profile_score_internal_3x1_nucleotides.begin(), profile_score_internal_3x1_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
    FillCounts(profile_score_internal_2x3_nucleotides.begin(), profile_score_internal_2x3_nucleotides.end(), RealT(0));
    FillCounts(profile_score_internal_3x2_nucleotides.begin(), profile_score_internal_3x2_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
    FillCounts(profile_score_internal_3x3_nucleotides.begin(), profile_score_internal_3x3_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
    FillCounts(profile_score_internal_1x4_nucleotides.begin(), profile_score_internal_1x4_nucleotides.end(), RealT(0));
    FillCounts(profile_score_internal_4x1_nucleotides.begin(), profile_score_internal_4x1_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
    FillCounts(profile_score_internal_2x4_nucleotides.begin(), profile_score_internal_2x4_nucleotides.end(), RealT(0));
    FillCounts(profile_score_internal_4x2_nucleotides.begin(), profile_score_internal_4x2_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
    FillCounts(profile_score_internal_3x4_nucleotides.begin(), profile_score_internal_3x4_nucleotides.end(), RealT(0));
    FillCounts(profile_score_internal_4x3_nucleotides.begin(), profile_score_internal_4x3_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
    FillCounts(profile_score_internal_4x4_nucleotides.begin(), profile_score_internal_4x4_nucleotides.end(), RealT(0));
#endif
    //自作終わり
#if PARAMS_HELIX_STACKING
    FillCounts(profile_score_helix_stacking.begin(), profile_score_helix_stacking.end(), RealT(0));
#endif
#if PARAMS_HELIX_CLOSING
    FillCounts(profile_score_helix_closing.begin(), profile_score_helix_closing.end(), RealT(0));
#endif
#if PARAMS_DANGLE
    FillCounts(profile_score_dangle_left.begin(), profile_score_dangle_left.end(), RealT(0));
    FillCounts(profile_score_dangle_right.begin(), profile_score_dangle_right.end(), RealT(0));
#endif

#endif

}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::FinalizeCounts()
//
// Apply any needed transformations to counts.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::FinalizeCounts()
{
    auto& pc = *parameter_count;
#if FAST_HELIX_LENGTHS

    // reverse helix partial sums    
    std::vector<std::pair<RealT,uint> > reverse_sums(cache_score_helix_sums);
    
    for (int i = 1; i <= L; i++)
    {
        for (int j = L; j >= i+3; j--)
        {
            // the "if" conditions here can be omitted
            
            if (allow_paired[offset[i+1]+j-1])
            {
                CountBasePair(i+1,j-1,reverse_sums[(i+j)*L+j-i].second);
                if (allow_paired[offset[i]+j])
                {
                    CountHelixStacking(i,j,reverse_sums[(i+j)*L+j-i].second);
                }
                else
                {
                    Assert(Abs(double(reverse_sums[(i+j)*L+j-i].second)) < 1e-8, "Should be zero.");
                }
            }
            else
            {
                Assert(Abs(double(reverse_sums[(i+j)*L+j-i-2].second)) < 1e-8, "Should be zero.");
            }
            
            reverse_sums[(i+j)*L+j-i-2].second += reverse_sums[(i+j)*L+j-i].second;
        }
    }
#endif

    // perform transformations
#if PARAMS_BASE_PAIR_DIST
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
        for (int j = BP_DIST_THRESHOLDS[i]; j <= BP_DIST_LAST_THRESHOLD; j++)
            pc.base_pair_dist_at_least(i) += cache_score_base_pair_dist[j].second;
#endif
    
#if PARAMS_HAIRPIN_LENGTH
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
        for (int j = i; j <= D_MAX_HAIRPIN_LENGTH; j++)
            pc.hairpin_length_at_least(i) += cache_score_hairpin_length[j].second;
#endif
    
#if PARAMS_HELIX_LENGTH
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
        for (int j = i; j <= D_MAX_HELIX_LENGTH; j++)
            pc.helix_length_at_least(i) += cache_score_helix_length[j].second;
#endif

    // allocate temporary storage
#if PARAMS_BULGE_LENGTH
    uint temp_cache_counts_bulge_length[D_MAX_BULGE_LENGTH+1];
    std::fill(temp_cache_counts_bulge_length, temp_cache_counts_bulge_length + D_MAX_BULGE_LENGTH+1, 0);
#endif
    
#if PARAMS_INTERNAL_LENGTH
    uint temp_cache_counts_internal_length[D_MAX_INTERNAL_LENGTH+1];
    std::fill(temp_cache_counts_internal_length, temp_cache_counts_internal_length + D_MAX_INTERNAL_LENGTH+1, 0);
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    uint temp_cache_counts_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
    std::fill(temp_cache_counts_internal_symmetric_length, temp_cache_counts_internal_symmetric_length + D_MAX_INTERNAL_SYMMETRIC_LENGTH+1, 0);
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    uint temp_cache_counts_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
    std::fill(temp_cache_counts_internal_asymmetry, temp_cache_counts_internal_asymmetry + D_MAX_INTERNAL_ASYMMETRY+1, 0);
#endif

    // compute contributions
    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
    {
        for (int l2 = 0; l1+l2 <= C_MAX_SINGLE_LENGTH; l2++)
        {
            // skip over stacking pairs
            if (l1 == 0 && l2 == 0) continue;

            // consider bulge loops
            if (l1 == 0 || l2 == 0)
            {
#if PARAMS_BULGE_LENGTH
                temp_cache_counts_bulge_length[std::min(D_MAX_BULGE_LENGTH, l1+l2)] += cache_score_single[l1][l2].second;
#endif
            }

            // consider internal loops
            else
            {
#if PARAMS_INTERNAL_EXPLICIT
                if (l1 <= D_MAX_INTERNAL_EXPLICIT_LENGTH && l2 <= D_MAX_INTERNAL_EXPLICIT_LENGTH)
                    pc.internal_explicit(l1, l2) += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_LENGTH
                temp_cache_counts_internal_length[std::min(D_MAX_INTERNAL_LENGTH, l1+l2)] += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_SYMMETRY
                if (l1 == l2)
                    temp_cache_counts_internal_symmetric_length[std::min(D_MAX_INTERNAL_SYMMETRIC_LENGTH, l1)] += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_ASYMMETRY
                temp_cache_counts_internal_asymmetry[std::min(D_MAX_INTERNAL_ASYMMETRY, Abs(l1-l2))] += cache_score_single[l1][l2].second;
#endif
            }
        }
    }

#if PARAMS_BULGE_LENGTH
    for (int i = 0; i <= D_MAX_BULGE_LENGTH; i++)
        for (int j = i; j <= D_MAX_BULGE_LENGTH; j++)
            pc.bulge_length_at_least(i) += temp_cache_counts_bulge_length[j];
#endif
    
#if PARAMS_INTERNAL_LENGTH
    for (int i = 0; i <= D_MAX_INTERNAL_LENGTH; i++)
        for (int j = i; j <= D_MAX_INTERNAL_LENGTH; j++)
            pc.internal_length_at_least(i) += temp_cache_counts_internal_length[j];
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    for (int i = 0; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
        for (int j = i; j <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; j++)
            pc.internal_symmetric_length_at_least(i) += temp_cache_counts_internal_symmetric_length[j];
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    for (int i = 0; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
        for (int j = i; j <= D_MAX_INTERNAL_ASYMMETRY; j++)
            pc.internal_asymmetry_at_least(i) += temp_cache_counts_internal_asymmetry[j];
#endif

    // finalize profile counts
#if PROFILE
    for (int i = 0; i <= L; i++)
    {
        for (int j = 0; j <= L; j++)
        {
#if PARAMS_BASE_PAIR
            {
                const int pos[2] = {i, j};
                ConvertProfileCount(profile_score_base_pair[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_base_pair));
            }
#endif
#if PARAMS_TERMINAL_MISMATCH
            {
                const int pos[4] = {i, j+1, i+1, j};
                ConvertProfileCount(profile_score_terminal_mismatch[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_terminal_mismatch));
            }
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ConvertProfileCount(profile_score_hairpin_3_nucleotides[i].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_3_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[4] = {i+1, i+2, i+3, i+4};
                ConvertProfileCount(profile_score_hairpin_4_nucleotides[i].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_4_nucleotides));
            }
#endif
	    //自作パラメータ
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[5] = {i+1, i+2, i+3, i+4,i+5};
                ConvertProfileCount(profile_score_hairpin_5_nucleotides[i].second, pos, 5, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_5_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[6] = {i+1, i+2, i+3, i+4, ,i+5, i+6};
                ConvertProfileCount(profile_score_hairpin_6_nucleotides[i].second, pos, 6, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_6_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[7] = {i+1, i+2, i+3, i+4, i+5, i+6, i+7};
                ConvertProfileCount(profile_score_hairpin_7_nucleotides[i].second, pos, 7, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_7_nucleotides));
            }
#endif
	    //自作終了
#if PARAMS_BULGE_0x1_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[1] = {j};
                ConvertProfileCount(profile_score_bulge_0x1_nucleotides[j].second, pos, 1, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x1_nucleotides));
            }
            if (j == 0)
            {
                const int pos[1] = {i+1};
                ConvertProfileCount(profile_score_bulge_1x0_nucleotides[i].second, pos, 1, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_1x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[2] = {j-1, j};
                ConvertProfileCount(profile_score_bulge_0x2_nucleotides[j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x2_nucleotides));
            }
            if (j == 0)
            {
                const int pos[2] = {i+1, i+2};
                ConvertProfileCount(profile_score_bulge_2x0_nucleotides[i].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_2x0_nucleotides));
            }
#endif            
#if PARAMS_BULGE_0x3_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[3] = {j-2, j-1, j};
                ConvertProfileCount(profile_score_bulge_0x3_nucleotides[j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x3_nucleotides));
            }
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ConvertProfileCount(profile_score_bulge_3x0_nucleotides[i].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_3x0_nucleotides));
            }
#endif
	    //自作パラメータ
#if PARAMS_BULGE_0x4_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[4] = {j-3,j-2, j-1, j};
                ConvertProfileCount(profile_score_bulge_0x4_nucleotides[j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x4_nucleotides));
            }
            if (j == 0)
            {
                const int pos[4] = {i+1, i+2, i+3, i+4};
                ConvertProfileCount(profile_score_bulge_4x0_nucleotides[i].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_4x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[5] = {j-4, j-3, j-2, j-1, j};
                ConvertProfileCount(profile_score_bulge_0x5_nucleotides[j].second, pos, 5, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x5_nucleotides));
            }
            if (j == 0)
            {
                const int pos[5] = {i+1, i+2, i+3,i+4,i+5};
                ConvertProfileCount(profile_score_bulge_5x0_nucleotides[i].second, pos, 5, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_5x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[6] = {j-5, j-4, j-3, j-2, j-1, j};
                ConvertProfileCount(profile_score_bulge_0x6_nucleotides[j].second, pos, 6, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x6_nucleotides));
            }
            if (j == 0)
            {
                const int pos[6] = {i+1, i+2, i+3, i+4, i+5, i+6};
                ConvertProfileCount(profile_score_bulge_6x0_nucleotides[i].second, pos, 6, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_6x0_nucleotides));
            }
#endif            
	    //自作終了
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
            {
                const int pos[2] = {i+1, j};
                ConvertProfileCount(profile_score_internal_1x1_nucleotides[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_1x1_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
            {
                const int pos[3] = {i+1, j-1, j};
                ConvertProfileCount(profile_score_internal_1x2_nucleotides[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_1x2_nucleotides));
            }
            {
                const int pos[3] = {i+1, i+2, j};
                ConvertProfileCount(profile_score_internal_2x1_nucleotides[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_2x1_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
            {
                const int pos[4] = {i+1, i+2, j-1, j};
                ConvertProfileCount(profile_score_internal_2x2_nucleotides[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_2x2_nucleotides));
            }
#endif     
	    //自作パラメータ
#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
            {
                const int pos[4] = {i+1, j-2, j-1, j};
                ConvertProfileCount(profile_score_internal_1x3_nucleotides[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_1x3_nucleotides));
            }
            {
                const int pos[4] = {i+1, i+2, i+3, j};
                ConvertProfileCount(profile_score_internal_3x1_nucleotides[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_3x1_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
            {
                const int pos[5] = {i+1, i+2, j-2, j-1, j};
                ConvertProfileCount(profile_score_internal_2x3_nucleotides[i*(L+1)+j].second, pos, 5, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_2x3_nucleotides));
            }
            {
                const int pos[5] = {i+1, i+2, i+3, j-1, j};
                ConvertProfileCount(profile_score_internal_3x2_nucleotides[i*(L+1)+j].second, pos, 5, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_3x2_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
            {
                const int pos[6] = {i+1, i+2, i+3, j-2, j-1, j};
                ConvertProfileCount(profile_score_internal_3x3_nucleotides[i*(L+1)+j].second, pos, 6, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_3x3_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
	    {
                const int pos[5] = {i+1, j-3, j-2, j-1, j};
                ConvertProfileCount(profile_score_internal_1x4_nucleotides[i*(L+1)+j].second, pos, 5, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_1x4_nucleotides));
	    }
	    {
                const int pos[5] = {i+1, i+2, i+3, i+4, j};
                ConvertProfileCount(profile_score_internal_4x1_nucleotides[i*(L+1)+j].second, pos, 5, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_4x1_nucleotides));
	    }
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
	    {
                const int pos[6] = {i+1, i+2, j-3, j-2, j-1, j};
                ConvertProfileCount(profile_score_internal_2x4_nucleotides[i*(L+1)+j].second, pos, 6, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_2x4_nucleotides));
	    }
	    {
                const int pos[6] = {i+1, i+2, i+3, i+4, j-1, j};
                ConvertProfileCount(profile_score_internal_4x2_nucleotides[i*(L+1)+j].second, pos, 6, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_4x2_nucleotides));
	    }
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
	    {
                const int pos[7] = {i+1, i+2, i+3, j-3, j-2, j-1, j};
                ConvertProfileCount(profile_score_internal_3x4_nucleotides[i*(L+1)+j].second, pos, 7, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_3x4_nucleotides));
	    }
	    {
                const int pos[7] = {i+1, i+2, i+3, i+4, j-2, j-1, j};
                ConvertProfileCount(profile_score_internal_4x3_nucleotides[i*(L+1)+j].second, pos, 7, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_4x3_nucleotides));
	    }
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
	    {
                const int pos[8] = {i+1, i+2, i+3, i+4, j-3, j-2, j-1, j};
                ConvertProfileCount(profile_score_internal_4x4_nucleotides[i*(L+1)+j].second, pos, 8, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_4x4_nucleotides));
	    }
#endif
	    
	    //自作終了
#if PARAMS_HELIX_STACKING
            {
                const int pos[4] = {i, j, i+1, j-1};
                ConvertProfileCount(profile_score_helix_stacking[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_helix_stacking));
            }
#endif
#if PARAMS_HELIX_CLOSING
            {
                const int pos[2] = {i, j+1};
                ConvertProfileCount(profile_score_helix_closing[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_helix_closing));
            }
#endif
#if PARAMS_DANGLE
            {
                const int pos[3] = {i, j+1, i+1};
                ConvertProfileCount(profile_score_dangle_left[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_dangle_left));
            }
            {
                const int pos[3] = {i, j+1, j};
                ConvertProfileCount(profile_score_dangle_right[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_dangle_right));
            }
#endif
        }
    }
#endif

}

#if PROFILE

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ConvertProfileCount()
//
// Convert profile count for a single location.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ConvertProfileCount(const RealT &profile_score, const int *pos, int dimensions, std::pair<RealT,RealT> *table)
{
    // consider all sequences
    for (int k = 0; k < N; k++)
    {
        bool valid = true;
        int index = 0;
        int *seq = &A[k*(L+1)];
        
        // extract letters of the pattern for the current sequence
        for (int d = 0; valid && d < dimensions; d++)
        {
            if (pos[d] < 1 || pos[d] > L)
                valid = false;
            else
            {
                BYTE c = seq[pos[d]];
                if (c == BYTE(alphabet.size()))
                    valid = false;
                else
                    index = index * (M+1) + c;
            }
        }

        // add contribution of pattern to score
        if (valid) table[index].second += weights[k] * profile_score;
    }
}
            
#endif

//////////////////////////////////////////////////////////////////////
// InferenceEngine::UseLoss()
//
// Use per-position loss.  A loss is incurred if true_mapping[i] !=
// UNKNOWN && solution[i] != true_mapping[i].
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::UseLoss(const std::vector<int> &true_mapping, RealT example_loss)
{
    Assert(int(true_mapping.size()) == L+1, "Mapping of incorrect length!");
    cache_initialized = false;
    
    // compute number of pairings
    int num_pairings = 0;
    for (int i = 1; i <= L; i++)
        if (true_mapping[i] != SStruct::UNKNOWN && true_mapping[i] != SStruct::UNPAIRED)
            ++num_pairings;

    //RealT per_position_loss = example_loss+(1/ RealT(num_pairings));
    RealT per_position_loss = example_loss;
    //std::cout << per_position_loss << std::endl;
    // RealT per_position_loss = 2*example_loss;
    
    // compute the penalty for each position that we declare to be unpaired
    for (int i = 0; i <= L; i++)
    {
        loss_unpaired_position[i] =
            ((i == 0 || true_mapping[i] == SStruct::UNKNOWN || true_mapping[i] == SStruct::UNPAIRED) ? RealT(0) : per_position_loss);
    }

    // now, compute the penalty for declaring ranges of positions to be unpaired;
    // also, compute the penalty for matching positions s[i] and s[j].
    for (int i = 0; i <= L; i++)
    {
        loss_unpaired[offset[i]+i] = RealT(0);
        loss_paired[offset[i]+i] = RealT(NEG_INF);
        for (int j = i+1; j <= L; j++)
        {
            loss_unpaired[offset[i]+j] = 
                loss_unpaired[offset[i]+j-1] +
                loss_unpaired_position[j];
            loss_paired[offset[i]+j] = 
                ((i == 0 || true_mapping[i] == SStruct::UNKNOWN || true_mapping[i] == SStruct::UNPAIRED || true_mapping[i] == j) ? RealT(0) : per_position_loss) +
                ((i == 0 || true_mapping[j] == SStruct::UNKNOWN || true_mapping[j] == SStruct::UNPAIRED || true_mapping[j] == i) ? RealT(0) : per_position_loss);
        }
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::UseConstraints()
//
// Use known secondary structure mapping.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::UseConstraints(const std::vector<int> &true_mapping)
{
    Assert(int(true_mapping.size()) == L+1, "Supplied mapping of incorrect length!");
    cache_initialized = false;
    
    // determine whether we allow each position to be unpaired
    for (int i = 1; i <= L; i++)
    { 
        allow_unpaired_position[i] =
            (true_mapping[i] == SStruct::UNKNOWN || 
             true_mapping[i] == SStruct::UNPAIRED);
    }

    // determine whether we allow ranges of positions to be unpaired;
    // also determine which base-pairings we allow
    for (int i = 0; i <= L; i++)
    {
        allow_unpaired[offset[i]+i] = 1;
        allow_paired[offset[i]+i] = 0;
        for (int j = i+1; j <= L; j++)
        {
            allow_unpaired[offset[i]+j] = 
                allow_unpaired[offset[i]+j-1] && 
                allow_unpaired_position[j];
            allow_paired[offset[i]+j] =
                (i > 0 &&
                 (true_mapping[i] == SStruct::UNKNOWN || true_mapping[i] == SStruct::PAIRED || true_mapping[i] == j) &&
                 (true_mapping[j] == SStruct::UNKNOWN || true_mapping[j] == SStruct::PAIRED || true_mapping[j] == i) &&
                 (allow_noncomplementary || IsComplementary(i,j)));
        }
    }
}



// score for leaving s[i] unpaired

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreUnpairedPosition(int i) const
{
#if defined(HAMMING_LOSS)
    return loss_unpaired_position[i]+reactivity_unpaired_position[i];
#else
    return RealT(0)+reactivity_unpaired_position[i];
#endif
}

template<class RealT>
inline void InferenceEngine<RealT>::CountUnpairedPosition(int i, uint v)
{
}

// score for leaving s[i+1...j] unpaired

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreUnpaired(int i, int j) const
{
#if defined(HAMMING_LOSS)
    return loss_unpaired[offset[i]+j]+reactivity_unpaired[offset[i]+j];
#else
    return RealT(0)+reactivity_unpaired[offset[i]+j];
#endif
}

template<class RealT>
inline void InferenceEngine<RealT>::CountUnpaired(int i,int j, uint v)
{
}

// score for a base pair which is not part of any helix

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreIsolated() const
{
#if PARAMS_ISOLATED_BASE_PAIR
    const auto& pm = *parameter_manager;
    return pm.isolated_base_pair();
#else
    return RealT(0);
#endif
}

template<class RealT>
inline void InferenceEngine<RealT>::CountIsolated(uint v)
{
#if PARAMS_ISOLATED_BASE_PAIR
    auto& pc = *parameter_count;
    pc.isolated_base_pair() += (v);
#endif
}

// base score for a multi-branch loop

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreMultiBase() const
{
#if PARAMS_MULTI_LENGTH
    const auto& pm = *parameter_manager;
    return pm.multi_base();
#else
    return RealT(0);
#endif    
}
    
template<class RealT>
inline void InferenceEngine<RealT>::CountMultiBase(uint v)
{
#if PARAMS_MULTI_LENGTH
    auto& pc = *parameter_count;
    pc.multi_base() += (v);
#endif
}

// score for a base-pair adjacent to a multi-branch loop

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreMultiPaired() const
{
#if PARAMS_MULTI_LENGTH
    const auto& pm = *parameter_manager;
    return pm.multi_paired();
#else
    return RealT(0);
#endif
}

template<class RealT>
inline void InferenceEngine<RealT>::CountMultiPaired(uint v)
{
#if PARAMS_MULTI_LENGTH
    auto& pc = *parameter_count;
    pc.multi_paired() += (v);
#endif
}

// score for each unpaired position in a multi-branch loop

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreMultiUnpaired(int i) const
{
#if PARAMS_MULTI_LENGTH
    const auto& pm = *parameter_manager;
    return pm.multi_unpaired() + ScoreUnpairedPosition(i);
#else
    return ScoreUnpairedPosition(i);
#endif
}

template<class RealT>
inline void InferenceEngine<RealT>::CountMultiUnpaired(int i, uint v)
{
#if PARAMS_MULTI_LENGTH
    auto& pc = *parameter_count;
    pc.multi_unpaired() += v; 
    CountUnpairedPosition(i,v);
#else
    CountUnpairedPosition(i,v);
#endif
}

// score for each base-pair adjacent to an external loop

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreExternalPaired() const
{
#if PARAMS_EXTERNAL_LENGTH
    const auto& pm = *parameter_manager;
    return pm.external_paired();
#else
    return RealT(0);
#endif
}

template<class RealT>
inline void InferenceEngine<RealT>::CountExternalPaired(uint v)
{
#if PARAMS_EXTERNAL_LENGTH
    auto& pc = *parameter_count;
    pc.external_paired() += (v);
#endif
}

// score for each unpaired position in an external loop

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreExternalUnpaired(int i) const
{
#if PARAMS_EXTERNAL_LENGTH
    const auto& pm = *parameter_manager;
    return pm.external_unpaired() + ScoreUnpairedPosition(i);
#else
    return ScoreUnpairedPosition(i);
#endif
}

template<class RealT>
inline void InferenceEngine<RealT>::CountExternalUnpaired(int i, uint v)
{
#if PARAMS_EXTERNAL_LENGTH
    auto& pc = *parameter_count;
    pc.external_unpaired() += v;
    CountUnpairedPosition(i,v);
#else
    CountUnpairedPosition(i,v); 
#endif
}

// score for a helix stacking pair of the form:
//
//       |         |
//    s[i+1] == s[j-1]
//       |         |
//     s[i] ==== s[j]
//       |         |

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreHelixStacking(int i, int j) const
{
#if PARAMS_HELIX_STACKING
#if PROFILE
    return profile_score_helix_stacking[i*(L+1)+j].first;
#else
    const auto& pm = *parameter_manager;
    return pm.helix_stacking(s[i], s[j], s[i+1], s[j-1]);
#endif
#else
    return RealT(0);
#endif
}

template<class RealT>
inline void InferenceEngine<RealT>::CountHelixStacking(int i,int j, uint v)
{
#if PARAMS_HELIX_STACKING
#if PROFILE
    profile_score_helix_stacking[i*(L+1)+j].second += (v);
#else
    auto& pc = *parameter_count;
    pc.helix_stacking(s[i], s[j], s[i+1], s[j-1]) += (v);
#endif
#endif
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreJunctionA()
// InferenceEngine::CountJunctionA()
//
// Returns the score for an asymmetric junction at positions i
// and j such that (i,j+1) are base-paired and (i+1,j) are free.
//
// In an RNA structure, this would look like
//
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
// Note that the difference between ScoreJunctionA() and
// ScoreJunctionB() is that the former applies to multi-branch
// loops whereas the latter is used for hairpin loops and
// single-branch loops.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreJunctionA(int i, int j) const
{
    // i and j must be bounded away from the edges so that s[i] and s[j+1]
    // refer to actual nucleotides.  To allow us to use this macro when
    // scoring the asymmetric junction for an exterior loop whose closing
    // base pair include the first and last nucleotides of the sequence,
    // we allow i to be as large as L and j to be as small as 0.
    
    Assert(0 < i && i <= L && 0 <= j && j < L, "Invalid indices.");
    const auto& pm = *parameter_manager;

    return
        RealT(0)
#if PARAMS_HELIX_CLOSING
#if PROFILE
        + profile_score_helix_closing[i*(L+1)+j].first
#else                                          
        + pm.helix_closing(s[i], s[j+1])
#endif
#endif
#if PARAMS_DANGLE
#if PROFILE
        + (i < L ? profile_score_dangle_left[i*(L+1)+j].first : RealT(0))
        + (j > 0 ? profile_score_dangle_right[i*(L+1)+j].first : RealT(0))
#else
        + (i < L ? pm.dangle_left(s[i], s[j+1], s[i+1]) : RealT(0))
        + (j > 0 ? pm.dangle_right(s[i], s[j+1], s[j]) : RealT(0))
#endif
#endif
        ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountJunctionA(int i, int j, uint value)
{
    Assert(0 < i && i <= L && 0 <= j && j < L, "Invalid indices.");
    auto& pc = *parameter_count;
    
#if PARAMS_HELIX_CLOSING
#if PROFILE
    profile_score_helix_closing[i*(L+1)+j].second += value;
#else
    pc.helix_closing(s[i], s[j+1]) += value;
#endif
#endif
#if PARAMS_DANGLE
#if PROFILE
    if (i < L) profile_score_dangle_left[i*(L+1)+j].second += value;
    if (j > 0) profile_score_dangle_right[i*(L+1)+j].second += value;
#else                                                               
    if (i < L) pc.dangle_left(s[i], s[j+1], s[i+1]) += value;
    if (j > 0) pc.dangle_right(s[i], s[j+1], s[j]) += value;
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreJunctionB()
// InferenceEngine::CountJunctionB()
//
// Returns the score for a symmetric junction at positions i
// and j such that (i,j+1) are base-paired and (i+1,j) are free.
//
// In an RNA structure, this would look like
//
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
// Note that the difference between ScoreJunctionA() and
// ScoreJunctionB() is that the former applies to multi-branch
// loops whereas the latter is used for hairpin loops and
// single-branch loops.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreJunctionB(int i, int j) const
{
    // The bounds here are similar to the asymmetric junction case, with
    // the main difference being that symmetric junctions are not allowed
    // for the exterior loop.  For this reason, i and j are bounded away
    // from the edges of the sequence (i.e., i < L && j > 0).
    //外部ループアリ
    Assert(0 < i && i < L && 0 < j && j < L, "Invalid indices.");
    const auto& pm = *parameter_manager;
    
    return RealT(0)
#if PARAMS_HELIX_CLOSING
#if PROFILE
        + profile_score_helix_closing[i*(L+1)+j].first
#else
        + pm.helix_closing(s[i], s[j+1])
#endif
#endif
#if PARAMS_TERMINAL_MISMATCH
#if PROFILE
        + profile_score_terminal_mismatch[i*(L+1)+j].first
#else                                           
        + pm.terminal_mismatch(s[i], s[j+1], s[i+1], s[j])
#endif
#endif
        ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountJunctionB(int i, int j, uint value)
{
    Assert(0 < i && i < L && 0 < j && j < L, "Invalid indices.");
    auto& pc = *parameter_count;
    
#if PARAMS_HELIX_CLOSING
#if PROFILE
    profile_score_helix_closing[i*(L+1)+j].second += value;
#else
    pc.helix_closing(s[i], s[j+1]) += value;
#endif
#endif
#if PARAMS_TERMINAL_MISMATCH
#if PROFILE
    profile_score_terminal_mismatch[i*(L+1)+j].second += value;
#else
    pc.terminal_mismatch(s[i], s[j+1], s[i+1], s[j]) += value;
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreBasePair()
// InferenceEngine::CountBasePair()
//
// Returns the score for a base-pairing between letters i and j.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreBasePair(int i, int j) const
{
    // Clearly, i and j must refer to actual letters of the sequence,
    // and no letter may base-pair to itself.
    
    Assert(0 < i && i <= L && 0 < j && j <= L && i != j, "Invalid base-pair");
    const auto& pm = *parameter_manager;
    
    return RealT(0)
#if defined(HAMMING_LOSS)
        + loss_paired[offset[i]+j]
#endif
#if PARAMS_BASE_PAIR
#if PROFILE
        + profile_score_base_pair[i*(L+1)+j].first
#else
        + pm.base_pair(s[i], s[j])
#endif
#endif
#if PARAMS_BASE_PAIR_DIST
        + cache_score_base_pair_dist[std::min(Abs(j - i), BP_DIST_LAST_THRESHOLD)].first
#endif
//pairのreactivity  
//#if DISCRETE==0
//                                  + use_reactivity*( sstruct.reactivity_pair[i] + sstruct.reactivity_pair[j])
//#endif
                                  ;
   // std::cout <<"reactivity in InferenceEngine.ipp is " << use_reactivity << std::endl; 
}

template<class RealT>
inline void InferenceEngine<RealT>::CountBasePair(int i, int j, uint value)
{
    Assert(0 < i && i <= L && 0 < j && j <= L && i != j, "Invalid base-pair");
    auto& pc = *parameter_count;
    
#if PARAMS_BASE_PAIR
#if PROFILE
    profile_score_base_pair[i*(L+1)+j].second += value;
#else
    pc.base_pair(s[i], s[j]) += value;
#endif
#endif
#if PARAMS_BASE_PAIR_DIST
    cache_score_base_pair_dist[std::min(Abs(j - i), BP_DIST_LAST_THRESHOLD)].second += value;
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreHairpin()
// InferenceEngine::CountHairpin()
//
// Returns the score for a hairpin spanning positions i to j.
//
// In an RNA structure, this would look like
//
//                           ...
//                       /         \. 
//                   x[i+2]       x[j-1]
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreHairpin(int i, int j) const
{
    // The constraints i > 0 && j < L ensure that s[i] and s[j+1] refer to
    // nucleotides which could base-pair.  The remaining constraint ensures
    // that only valid hairpins are considered.
    
    Assert(0 < i && i + C_MIN_HAIRPIN_LENGTH <= j && j < L, "Hairpin boundaries invalid.");
    const auto& pm = *parameter_manager;
    
    return 
        ScoreUnpaired(i,j)
        + ScoreJunctionB(i,j)
#if PARAMS_HAIRPIN_LENGTH
        + cache_score_hairpin_length[std::min(j - i, D_MAX_HAIRPIN_LENGTH)].first
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
#if PROFILE
        + (j - i == 3 ? profile_score_hairpin_3_nucleotides[i].first : RealT(0))
#else
        + (j - i == 3 ? pm.hairpin_3_nucleotides(s[i+1], s[i+2], s[i+3]) : RealT(0))
#endif                                          
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
#if PROFILE
        + (j - i == 4 ? profile_score_hairpin_4_nucleotides[i].first : RealT(0))
#else
        + (j - i == 4 ? pm.hairpin_4_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4]) : RealT(0))
#endif
#endif
      //自作パラメータ
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
#if PROFILE
      + (j - i == 5 ? profile_score_hairpin_5_nucleotides[i].first : RealT(0))
#else
      + (j - i == 5 ? pm.hairpin_5_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5]) : RealT(0))
#endif
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
#if PROFILE
      + (j - i == 6 ? profile_score_hairpin_6_nucleotides[i].first : RealT(0))
#else
      + (j - i == 6 ? pm.hairpin_6_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5], s[i+6]) : RealT(0))
#endif
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
#if PROFILE
      + (j - i == 7 ? profile_score_hairpin_7_nucleotides[i].first : RealT(0))
#else
      + (j - i == 7 ? pm.hairpin_7_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5], s[i+6], s[i+7]) : RealT(0))
#endif
#endif
      //自作終了
      ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountHairpin(int i, int j, uint value)
{
    Assert(0 < i && i + C_MIN_HAIRPIN_LENGTH <= j && j < L, "Hairpin boundaries invalid.");
    auto& pc = *parameter_count;
    
    CountUnpaired(i,j,value);
    CountJunctionB(i,j,value);
#if PARAMS_HAIRPIN_LENGTH
    cache_score_hairpin_length[std::min(j - i, D_MAX_HAIRPIN_LENGTH)].second += value;
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
#if PROFILE
    if (j - i == 3) profile_score_hairpin_3_nucleotides[i].second += value;
#else
    if (j - i == 3) pc.hairpin_3_nucleotides(s[i+1], s[i+2], s[i+3]) += value;
#endif
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
#if PROFILE
    if (j - i == 4) profile_score_hairpin_4_nucleotides[i].second += value;
#else
    if (j - i == 4) pc.hairpin_4_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4]) += value;
#endif
#endif
    //自作パラメータ
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
#if PROFILE
    if (j - i == 5) profile_score_hairpin_5_nucleotides[i].second += value;
#else
    if (j - i == 5) pc.hairpin_5_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5]) += value;
#endif
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
#if PROFILE
    if (j - i == 6) profile_score_hairpin_6_nucleotides[i].second += value;
#else
    if (j - i == 6) pc.hairpin_6_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5], s[i+6]) += value;
#endif
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
#if PROFILE
    if (j - i == 7) profile_score_hairpin_7_nucleotides[i].second += value;
#else
    if (j - i == 7) score_hairpin_7_nucleotides[s[i+1]][s[i+2]][s[i+3]][s[i+4]][s[i+5]][s[i+6]][s[i+7]].second += value;
#endif
#endif
    //自作終了
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreHelix()
// InferenceEngine::CountHelix()
//
// Returns the score for a helix of length m starting at positions
// i and j.  All base-pairs except for x[i+1]-x[j] are scored.
//
// In an RNA structure, this would look like
//
//                           ...
//                       \          /
// position i+m ------->  o        o  <----- position j-m
//                     x[i+3] -- x[j-2]
//                        |        |
//                     x[i+2] -- x[j-1]
//                        |        |
//                     x[i+1] -- x[j]
// position i --------->  o        o  <----- position j
//                       /          \.
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreHelix(int i, int j, int m) const
{
    // First, i >= 0 && j <= L are obvious sanity-checks to make sure that
    // things are within range.  The check that i+2*m <= j ensures that there
    // are enough nucleotides to allow a helix of length m.
    
    Assert(0 <= i && i + 2 * m <= j && j <= L, "Helix boundaries invalid.");
    Assert(2 <= m && m <= D_MAX_HELIX_LENGTH, "Helix length invalid.");
    
#if FAST_HELIX_LENGTHS
    
    return
        cache_score_helix_sums[(i+j+1)*L+j-i-1].first - cache_score_helix_sums[(i+j+1)*L+j-i-m-m+1].first
#if PARAMS_HELIX_LENGTH
        + cache_score_helix_length[m].first
#endif
        ;
    
#else 
    
    RealT ret = RealT(0);
    for (int k = 1; k < m; k++)
        ret += ScoreHelixStacking(i+k,j-k+1) + ScoreBasePair(i+k+1,j-k);
    
#if PARAMS_HELIX_LENGTH
    ret += cache_score_helix_length[m].first;
#endif

    return ret;

#endif
  
}

template<class RealT>
inline void InferenceEngine<RealT>::CountHelix(int i, int j, int m, uint value)
{
    Assert(0 <= i && i + 2 * m <= j && j <= L, "Helix boundaries invalid.");
    Assert(2 <= m && m <= D_MAX_HELIX_LENGTH, "Helix length invalid.");
    
#if FAST_HELIX_LENGTHS
    
    cache_score_helix_sums[(i+j+1)*L+j-i-1].second += value;
    cache_score_helix_sums[(i+j+1)*L+j-i-m-m+1].second -= value;
    
#else
    
    for (int k = 1; k < m; k++)
    {
        CountHelixStacking(i+k,j-k+1,value);
        CountBasePair(i+k+1,j-k,value);
    }
    
#endif
    
#if PARAMS_HELIX_LENGTH
    cache_score_helix_length[m].second += value;
#endif
    
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreSingleNucleotides()
// InferenceEngine::CountSingleNucleotides()
//
// Returns the score for nucleotides in a single-branch loop 
// spanning i to j and p to q.
//
// In an RNA structure, this would look like
//
//                       ...      ...
//                        |        |
//                      x[p+1] -- x[q]
// position p -------->  o          o  <----- position q
//                    x[p]        x[q+1]
//                      |            |
//                     ...          ...
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreSingleNucleotides(int i, int j, int p, int q) const
{
    // Nucleotides s[i] and s[j+1] must exist, hence the conditions i > 0 and j < L.
    // the condition p+2 <= q comes from the fact that there must be enough room for
    // at least one nucleotide on the other side of the single-branch loop.  This
    // loop should only be used for dealing with single-branch loops, not stacking pairs.
    
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");
    const auto& pm = *parameter_manager;

#if (!defined(NDEBUG) || PARAMS_BULGE_0x1_NUCLEOTIDES || PARAMS_BULGE_0x2_NUCLEOTIDES || PARAMS_BULGE_0x3_NUCLEOTIDES ||  PARAMS_BULGE_0x4_NUCLEOTIDES ||  PARAMS_BULGE_0x5_NUCLEOTIDES ||  PARAMS_BULGE_0x6_NUCLEOTIDES ||PARAMS_INTERNAL_1x1_NUCLEOTIDES || PARAMS_INTERNAL_1x2_NUCLEOTIDES || PARAMS_INTERNAL_2x2_NUCLEOTIDES || PARAMS_INTERNAL_1x3_NUCLEOTIDES || PARAMS_INTERNAL_2x3_NUCLEOTIDES || PARAMS_INTERNAL_3x3_NUCLEOTIDES || PARAMS_INTERNAL_1x4_NUCLEOTIDES || PARAMS_INTERNAL_2x4_NUCLEOTIDES || PARAMS_INTERNAL_3x4_NUCLEOTIDES || PARAMS_INTERNAL_4x4_NUCLEOTIDES )
    const int l1 = p - i;
    const int l2 = j - q;
    
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
#endif
    
    return 
        ScoreUnpaired(i,p)
        + ScoreUnpaired(q,j)
#if PARAMS_BULGE_0x1_NUCLEOTIDES
#if PROFILE
        + (l1 == 0 && l2 == 1 ? profile_score_bulge_0x1_nucleotides[j].first : RealT(0))
        + (l1 == 1 && l2 == 0 ? profile_score_bulge_1x0_nucleotides[i].first : RealT(0))
#else
        + (l1 == 0 && l2 == 1 ? pm.bulge_0x1_nucleotides(s[j]) : RealT(0))
        + (l1 == 1 && l2 == 0 ? pm.bulge_1x0_nucleotides(s[i+1]) : RealT(0))
#endif
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
#if PROFILE
        + (l1 == 0 && l2 == 2 ? profile_score_bulge_0x2_nucleotides[j].first : RealT(0))
        + (l1 == 2 && l2 == 0 ? profile_score_bulge_2x0_nucleotides[i].first : RealT(0))
#else
        + (l1 == 0 && l2 == 2 ? pm.bulge_0x2_nucleotides(s[j-1], s[j]) : RealT(0))
        + (l1 == 2 && l2 == 0 ? pm.bulge_2x0_nucleotides(s[i+1], s[i+2]): RealT(0))
#endif
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
#if PROFILE
        + (l1 == 0 && l2 == 3 ? profile_score_bulge_0x3_nucleotides[j].first : RealT(0))
        + (l1 == 3 && l2 == 0 ? profile_score_bulge_3x0_nucleotides[i].first : RealT(0))
#else
        + (l1 == 0 && l2 == 3 ? pm.bulge_0x3_nucleotides(s[j-2], s[j-1], s[j]) : RealT(0))
        + (l1 == 3 && l2 == 0 ? pm.bulge_3x0_nucleotides(s[i+1], s[i+2], s[i+3]) : RealT(0))
#endif
#endif
      //自作パラメータ
#if PARAMS_BULGE_0x4_NUCLEOTIDES
#if PROFILE
      + (l1 == 0 && l2 == 4 ? profile_score_bulge_0x4_nucleotides[j].first : RealT(0))
      + (l1 == 4 && l2 == 0 ? profile_score_bulge_4x0_nucleotides[i].first : RealT(0))
#else
      + (l1 == 0 && l2 == 4 ? pm.bulge_0x4_nucleotides(s[j-3], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 4 && l2 == 0 ? pm.bulge_4x0_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4]) : RealT(0))
#endif
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
#if PROFILE
      + (l1 == 0 && l2 == 5 ? profile_score_bulge_0x5_nucleotides[j].first : RealT(0))
      + (l1 == 5 && l2 == 0 ? profile_score_bulge_5x0_nucleotides[i].first : RealT(0))
#else
      + (l1 == 0 && l2 == 5 ? pm.bulge_0x5_nucleotides(s[j-4], s[j-3], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 5 && l2 == 0 ? pm.bulge_5x0_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5]) : RealT(0))
#endif
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
#if PROFILE
      + (l1 == 0 && l2 == 6 ? profile_score_bulge_0x6_nucleotides[j].first : RealT(0))
      + (l1 == 6 && l2 == 0 ? profile_score_bulge_6x0_nucleotides[i].first : RealT(0))
#else
      + (l1 == 0 && l2 == 6 ? pm.bulge_0x6_nucleotides(s[j-5], s[j-4], s[j-3], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 6 && l2 == 0 ? pm.bulge_6x0_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5], s[i+6]) : RealT(0))
#endif
#endif
#if PARAMS_BULGE_0x7_NUCLEOTIDES
#if PROFILE
      + (l1 == 0 && l2 == 7 ? profile_score_bulge_0x7_nucleotides[j].first : RealT(0))
      + (l1 == 7 && l2 == 0 ? profile_score_bulge_7x0_nucleotides[i].first : RealT(0))
#else
      + (l1 == 0 && l2 == 7 ? pm.bulge_0x7_nucleotides(s[j-6], s[j-5], s[j-4], s[j-3], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 7 && l2 == 0 ? pm.bulge_7x0_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5], s[i+6], s[i+7]) : RealT(0))
#endif
#endif
      //自作終了
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
#if PROFILE
        + (l1 == 1 && l2 == 1 ? profile_score_internal_1x1_nucleotides[i*(L+1)+j].first : RealT(0))
#else
        + (l1 == 1 && l2 == 1 ? pm.internal_1x1_nucleotides(s[i+1], s[j]) : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
#if PROFILE
        + (l1 == 1 && l2 == 2 ? profile_score_internal_1x2_nucleotides[i*(L+1)+j].first : RealT(0))
        + (l1 == 2 && l2 == 1 ? profile_score_internal_2x1_nucleotides[i*(L+1)+j].first : RealT(0))
#else
        + (l1 == 1 && l2 == 2 ? pm.internal_1x2_nucleotides(s[i+1], s[j-1], s[j]) : RealT(0))
        + (l1 == 2 && l2 == 1 ? pm.internal_2x1_nucleotides(s[i+1], s[i+2], s[j]) : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
#if PROFILE
        + (l1 == 2 && l2 == 2 ? profile_score_internal_2x2_nucleotides[i*(L+1)+j].first : RealT(0))
#else                                                                   
        + (l1 == 2 && l2 == 2 ? pm.internal_2x2_nucleotides(s[i+1], s[i+2], s[j-1], s[j]) : RealT(0))
#endif
#endif
      //自作パラメータ
#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
#if PROFILE
      + (l1 == 1 && l2 == 3 ? profile_score_internal_1x3_nucleotides[i*(L+1)+j].first : RealT(0))
      + (l1 == 3 && l2 == 1 ? profile_score_internal_3x1_nucleotides[i*(L+1)+j].first : RealT(0))
#else
      + (l1 == 1 && l2 == 3 ? pm.internal_1x3_nucleotides(s[i+1], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 3 && l2 == 1 ? pm.internal_3x1_nucleotides(s[i+1], s[i+2], s[i+3], s[j]) : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
#if PROFILE
      + (l1 == 2 && l2 == 3 ? profile_score_internal_2x3_nucleotides[i*(L+1)+j].first : RealT(0))
      + (l1 == 3 && l2 == 2 ? profile_score_internal_3x2_nucleotides[i*(L+1)+j].first : RealT(0))
#else
      + (l1 == 2 && l2 == 3 ? pm.internal_2x3_nucleotides(s[i+1], s[i+2], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 3 && l2 == 2 ? pm.internal_3x2_nucleotides(s[i+1], s[i+2], s[i+3], s[j-1], s[j]) : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
#if PROFILE
      + (l1 == 3 && l2 == 3 ? profile_score_internal_3x3_nucleotides[i*(L+1)+j].first : RealT(0))
#else
      + (l1 == 3 && l2 == 3 ? pm.internal_3x3_nucleotides(s[i+1], s[i+2], s[i+3], s[j-2], s[j-1], s[j]) : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
#if PROFILE
      + (l1 == 1 && l2 == 4 ? profile_score_internal_1x4_nucleotides[i*(L+1)+j].first : RealT(0))
      + (l1 == 4 && l2 == 1 ? profile_score_internal_4x1_nucleotides[i*(L+1)+j].first : RealT(0))
      #else
      + (l1 == 1 && l2 == 4 ? pm.internal_1x4_nucleotides(s[i+1], s[j-3], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 4 && l2 == 1 ? pm.internal_4x1_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[j]) : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
#if PROFILE
      + (l1 == 2 && l2 == 4 ? profile_score_internal_2x4_nucleotides[i*(L+1)+j].first : RealT(0))
      + (l1 == 4 && l2 == 2 ? profile_score_internal_4x2_nucleotides[i*(L+1)+j].first : RealT(0))
#else
      + (l1 == 2 && l2 == 4 ? pm.internal_2x4_nucleotides(s[i+1], s[i+2], s[j-3], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 4 && l2 == 2 ? pm.internal_4x2_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[j-1], s[j]) : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
#if PROFILE
      + (l1 == 3 && l2 == 4 ? profile_score_internal_3x4_nucleotides[i*(L+1)+j].first : RealT(0))
      + (l1 == 4 && l2 == 3 ? profile_score_internal_4x3_nucleotides[i*(L+1)+j].first : RealT(0))
#else
      + (l1 == 3 && l2 == 4 ? pm.internal_3x4_nucleotides(s[i+1], s[i+2], s[i+3], s[j-3], s[j-2], s[j-1], s[j]) : RealT(0))
      + (l1 == 4 && l2 == 3 ? pm.internal_4x3_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[j-2], s[j-1], s[j]) : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
#if PROFILE
      + (l1 == 4 && l2 == 4 ? profile_score_internal_4x4_nucleotides[i*(L+1)+j].first : RealT(0))
#else
      + (l1 == 4 && l2 == 4 ? pm.internal_4x4_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[j-3], s[j-2], s[j-1], s[j]) : RealT(0))
#endif
#endif      

      //自作終了
      ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountSingleNucleotides(int i, int j, int p, int q, uint value)
{
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");
    auto& pc = *parameter_count;
    
#if (!defined(NDEBUG) || PARAMS_BULGE_0x1_NUCLEOTIDES || PARAMS_BULGE_0x2_NUCLEOTIDES || PARAMS_BULGE_0x3_NUCLEOTIDES || PARAMS_BULGE_0x4_NUCLEOTIDES ||  PARAMS_BULGE_0x5_NUCLEOTIDES ||  PARAMS_BULGE_0x6_NUCLEOTIDES || PARAMS_INTERNAL_1x1_NUCLEOTIDES || PARAMS_INTERNAL_1x2_NUCLEOTIDES || PARAMS_INTERNAL_2x2_NUCLEOTIDES || PARAMS_INTERNAL_1x3_NUCLEOTIDES || PARAMS_INTERNAL_2x3_NUCLEOTIDES || PARAMS_INTERNAL_3x3_NUCLEOTIDES || PARAMS_INTERNAL_1x4_NUCLEOTIDES || PARAMS_INTERNAL_2x4_NUCLEOTIDES || PARAMS_INTERNAL_3x4_NUCLEOTIDES || PARAMS_INTERNAL_4x4_NUCLEOTIDES)
    const int l1 = p - i;
    const int l2 = j - q;
    
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
#endif
    
    CountUnpaired(i,p,value);
    CountUnpaired(q,j,value);
#if PARAMS_BULGE_0x1_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 1) profile_score_bulge_0x1_nucleotides[j].second += value;
    if (l1 == 1 && l2 == 0) profile_score_bulge_1x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 1) pc.bulge_0x1_nucleotides(s[j]) += value;
    if (l1 == 1 && l2 == 0) pc.bulge_1x0_nucleotides(s[i+1]) += value;
#endif
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 2) profile_score_bulge_0x2_nucleotides[j].second += value;
    if (l1 == 2 && l2 == 0) profile_score_bulge_2x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 2) pc.bulge_0x2_nucleotides(s[j-1], s[j]) += value;
    if (l1 == 2 && l2 == 0) pc.bulge_2x0_nucleotides(s[i+1], s[i+2]) += value;
#endif
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 3) profile_score_bulge_0x3_nucleotides[j].second += value;
    if (l1 == 3 && l2 == 0) profile_score_bulge_3x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 3) pc.bulge_0x3_nucleotides(s[j-2], s[j-1], s[j]) += value;
    if (l1 == 3 && l2 == 0) pc.bulge_3x0_nucleotides(s[i+1], s[i+2], s[i+3]) += value;
#endif
#endif
    //自作パラメータ
#if PARAMS_BULGE_0x4_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 4) profile_score_bulge_0x4_nucleotides[j].second += value;
    if (l1 == 4 && l2 == 0) profile_score_bulge_4x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 4) pc.bulge_0x4_nucleotides(s[j-3], s[j-2], s[j-1], s[j]) += value;
    if (l1 == 4 && l2 == 0) pc.bulge_4x0_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4]) += value;
#endif
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 5) profile_score_bulge_0x5_nucleotides[j].second += value;
    if (l1 == 5 && l2 == 0) profile_score_bulge_5x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 5) pc.bulge_0x5_nucleotides(s[j-4], s[j-3], s[j-2], s[j-1], s[j]) += value;
    if (l1 == 5 && l2 == 0) pc.bulge_5x0_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5]) += value;
#endif
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 6) profile_score_bulge_0x6_nucleotides[j].second += value;
    if (l1 == 6 && l2 == 0) profile_score_bulge_6x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 6) pc.bulge_0x6_nucleotides(s[j-5], s[j-4], s[j-3], s[j-2], s[j-1], s[j]) += value;
    if (l1 == 6 && l2 == 0) pc.bulge_6x0_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[i+5], s[i+6]) += value;
#endif
#endif
    //自作終了
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
#if PROFILE
    if (l1 == 1 && l2 == 1) profile_score_internal_1x1_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 1 && l2 == 1) pc.internal_1x1_nucleotides(s[i+1], s[j]) += value;
#endif
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
#if PROFILE
    if (l1 == 1 && l2 == 2) profile_score_internal_1x2_nucleotides[i*(L+1)+j].second += value;
    if (l1 == 2 && l2 == 1) profile_score_internal_2x1_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 1 && l2 == 2) pc.internal_1x2_nucleotides(s[i+1], s[j-1], s[j]) += value;
    if (l1 == 2 && l2 == 1) pc.internal_2x1_nucleotides(s[i+1], s[i+2], s[j]) += value;
#endif    
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
#if PROFILE
    if (l1 == 2 && l2 == 2) profile_score_internal_2x2_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 2 && l2 == 2) pc.internal_2x2_nucleotides(s[i+1], s[i+2], s[j-1], s[j]) += value;
#endif
#endif
    //自作パラメータ
#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
#if PROFILE
    if (l1 == 1 && l2 == 3) profile_score_internal_1x3_nucleotides[i*(L+1)+j].second += value;
    if (l1 == 3 && l2 == 1) profile_score_internal_3x1_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 1 && l2 == 3) pc.internal_1x3_nucleotides(s[i+1], s[j-2], s[j-1], s[j]) += value;
    if (l1 == 3 && l2 == 1) pc.internal_3x1_nucleotides(s[i+1], s[i+2], s[i+3], s[j]) += value;
#endif
#endif
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
#if PROFILE
    if (l1 == 2 && l2 == 3) profile_score_internal_2x3_nucleotides[i*(L+1)+j].second += value;
    if (l1 == 3 && l2 == 2) profile_score_internal_3x2_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 2 && l2 == 3) pc.internal_2x3_nucleotides(s[i+1], s[i+2], s[j-2], s[j-1], s[j]) += value;
    if (l1 == 3 && l2 == 2) pc.internal_3x2_nucleotides(s[i+1], s[i+2], s[i+3], s[j-1], s[j]) += value;
#endif
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
#if PROFILE
    if (l1 == 3 && l2 == 3) profile_score_internal_3x3_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 3 && l2 == 3) pc.internal_3x3_nucleotides(s[i+1], s[i+2], s[i+3], s[j-2], s[j-1], s[j]) += value;
#endif
#endif
#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
#if PROFILE
    if (l1 == 1 && l2 == 4) profile_score_internal_1x4_nucleotides[i*(L+1)+j].second += value;
    if (l1 == 4 && l2 == 1) profile_score_internal_4x1_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 1 && l2 == 4) pc.internal_1x4_nucleotides(s[i+1], s[j-3], s[j-2], s[j-1], s[j]) += value;
    if (l1 == 4 && l2 == 1) pc.internal_4x1_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[j]) += value;
#endif
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
#if PROFILE
    if (l1 == 2 && l2 == 4) profile_score_internal_2x4_nucleotides[i*(L+1)+j].second += value;
    if (l1 == 4 && l2 == 2) profile_score_internal_4x2_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 2 && l2 == 4) pc.internal_2x4_nucleotides(s[i+1], s[i+2], s[j-3], s[j-2], s[j-1], s[j]) += value;
    if (l1 == 4 && l2 == 2) pc.internal_4x2_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[j-1], s[j]) += value;
#endif
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
#if PROFILE
    if (l1 == 3 && l2 == 4) profile_score_internal_3x4_nucleotides[i*(L+1)+j].second += value;
    if (l1 == 4 && l2 == 3) profile_score_internal_4x3_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 3 && l2 == 4) pc.internal_3x4_nucleotides(s[i+1], s[i+2], s[i+3], s[j-3], s[j-2], s[j-1], s[j]) += value;
    if (l1 == 4 && l2 == 3) pc.internal_4x3_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[j-2], s[j-1], s[j]) += value;
#endif
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
#if PROFILE
    if (l1 == 4 && l2 == 4) profile_score_internal_4x4_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 4 && l2 == 4) pc.internal_4x4_nucleotides(s[i+1], s[i+2], s[i+3], s[i+4], s[j-3], s[j-2], s[j-1], s[j]) += value;
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreSingle()
// InferenceEngine::CountSingle()
//
// Returns the score for a single-branch loop spanning i to j and
// p to q.
//
// In an RNA structure, this would look like
//
//                       ...      ...
//                        |        |
//                      x[p+1] -- x[q]
// position p -------->  o          o  <----- position q
//                    x[p]        x[q+1]
//                      |            |
//                     ...          ...
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreSingle(int i, int j, int p, int q) const
{
    const int l1 = p - i;
    const int l2 = j - q;
    
    // Nucleotides s[i] and s[j+1] must exist, hence the conditions i > 0 and j < L.
    // the condition p+2 <= q comes from the fact that there must be enough room for
    // at least one nucleotide on the other side of the single-branch loop.  This
    // loop should only be used for dealing with single-branch loops, not stacking pairs.
    
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
    
    return 
        cache_score_single[l1][l2].first
        + ScoreBasePair(p+1,q)
        + ScoreJunctionB(i,j) 
        + ScoreJunctionB(q,p)
        + ScoreSingleNucleotides(i,j,p,q);
}

template<class RealT>
inline void InferenceEngine<RealT>::CountSingle(int i, int j, int p, int q, uint value)
{
    const int l1 = p - i;
    const int l2 = j - q;
    
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
    
    cache_score_single[l1][l2].second += value;
    CountBasePair(p+1,q,value);
    CountJunctionB(i,j,value);
    CountJunctionB(q,p,value);
    CountSingleNucleotides(i,j,p,q,value);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::EncodeTraceback()
// InferenceEngine::DecodeTraceback()
//
// Encode and decode traceback as an integer.  Here, i encodes
// a traceback type, and j encodes a length.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline int InferenceEngine<RealT>::EncodeTraceback(int i, int j) const
{
    Assert(0 <= i && i < NUM_TRACEBACK_TYPES && j >= 0, "Invalid values to encode as traceback.");
    return (j * NUM_TRACEBACK_TYPES) + i;
}

template<class RealT>
inline std::pair<int,int> InferenceEngine<RealT>::DecodeTraceback(int s) const
{
    return std::make_pair (s % NUM_TRACEBACK_TYPES, s / NUM_TRACEBACK_TYPES);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeViterbi()
//
// Run Viterbi algorithm.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputeViterbi()
{
    InitializeCache();
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif

#if CANDIDATE_LIST
    std::vector<int> candidates;
    candidates.reserve(L+1);
    long long int candidates_seen = 0;
    long long int candidates_possible = 0;
#endif
    
    // initialization

    F5t.clear(); F5t.resize(L+1, -1);
    FCt.clear(); FCt.resize(SIZE, -1);
    FMt.clear(); FMt.resize(SIZE, -1);
    FM1t.clear(); FM1t.resize(SIZE, -1);

    F5v.clear(); F5v.resize(L+1, RealT(NEG_INF));
    FCv.clear(); FCv.resize(SIZE, RealT(NEG_INF));
    FMv.clear(); FMv.resize(SIZE, RealT(NEG_INF));
    FM1v.clear(); FM1v.resize(SIZE, RealT(NEG_INF));
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    FEt.clear(); FEt.resize(SIZE, -1);
    FNt.clear(); FNt.resize(SIZE, -1);
    FEv.clear(); FEv.resize(SIZE, RealT(NEG_INF));
    FNv.clear(); FNv.resize(SIZE, RealT(NEG_INF));
#endif
    
    for (int i = L; i >= 0; i--)
    {
        
#if CANDIDATE_LIST
        candidates.clear();
#endif
        
        for (int j = i; j <= L; j++)
        {
            // FM2[i,j] = MAX (i<k<j : FM1[i,k] + FM[k,j])

            RealT FM2v = RealT(NEG_INF);
            int FM2t = -1;
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                UPDATE_MAX(FM2v, FM2t, FM1v[offset[i]+k] + FMv[offset[k]+j], k);
            
#else
            
#if !CANDIDATE_LIST
            
            if (i+2 <= j)
            {
                RealT *p1 = &(FM1v[offset[i]+i+1]);
                RealT *p2 = &(FMv[offset[i+1]+j]);
                for (int k = i+1; k < j; k++)
                {
                    UPDATE_MAX(FM2v, FM2t, (*p1) + (*p2), k);
                    ++p1;
                    p2 += L-k;
                }
            }
            
#else
            
            for (size_t kp = 0; kp < candidates.size(); kp++)
            {
                const int k = candidates[kp];
                UPDATE_MAX(FM2v, FM2t, FM1v[offset[i]+k] + FMv[offset[k]+j], k);
            }
            
            candidates_seen += (long long int) candidates.size();
            candidates_possible += (long long int) std::max(j-i-1,0);
            
#endif
            
#endif

#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
      
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = MAX [ScoreHairpin(i,j),
            //                MAX (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + MAX (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    UPDATE_MAX(best_v, best_t, ScoreHairpin(i,j), EncodeTraceback(TB_FN_HAIRPIN,0));
                
                // compute MAX (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        UPDATE_MAX(best_v, best_t,
                                   ScoreSingle(i,j,p,q) + FCv[offset[p+1]+q-1],
                                   EncodeTraceback(TB_FN_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q));
                    }
                }
                
#else
                
                {
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    int bp = -1, bq = -1;
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCv[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT score = (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) +
                                           ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                            if (score > best_v)
                            {
                                best_v = score;
                                bp = p;
                                bq = q;
                            }
                        }
                    }
                    
                    if (bp != -1 && bq != -1)
                        best_t = EncodeTraceback(TB_FN_SINGLE,(bp-i)*(C_MAX_SINGLE_LENGTH+1)+j-bq);
                }
#endif
                
                // compute MAX (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                UPDATE_MAX(best_v, best_t,
                           FM2v + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase(), 
                           EncodeTraceback(TB_FN_BIFURCATION,FM2t));
                
                FNv[offset[i]+j] = best_v;
                FNt[offset[i]+j] = best_t;
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = MAX [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    UPDATE_MAX(best_v, best_t, 
                               ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FEv[offset[i+1]+j-1],
                               EncodeTraceback(TB_FE_STACKING,0));
                }
                
                // compute FN(i,j)
                
                UPDATE_MAX(best_v, best_t, FNv[offset[i]+j], EncodeTraceback(TB_FE_FN,0));
                
                FEv[offset[i]+j] = best_v;
                FEt[offset[i]+j] = best_t;
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = MAX [ScoreIsolated() + FN(i,j),
            //                MAX (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute ScoreIsolated() + FN(i,j)
                
                UPDATE_MAX(best_v, best_t, ScoreIsolated() + FNv[offset[i]+j], EncodeTraceback(TB_FC_FN,0));
                
                // compute MAX (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    UPDATE_MAX(best_v, best_t, ScoreHelix(i-1,j+1,k) + FNv[offset[i+k-1]+j-k+1], EncodeTraceback(TB_FC_HELIX,k));
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        UPDATE_MAX(best_v, best_t, ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) +
                                   FEv[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1],
                                   EncodeTraceback(TB_FC_FE,0));
                }
                FCv[offset[i]+j] = best_v;
                FCt[offset[i]+j] = best_t;
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = MAX [ScoreHairpin(i,j),
            //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + MAX (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    UPDATE_MAX(best_v, best_t, ScoreHairpin(i,j), EncodeTraceback(TB_FC_HAIRPIN,0));
                
                // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])

#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        
                        UPDATE_MAX(best_v, best_t,
                                   FCv[offset[p+1]+q-1] +
                                   (p == i && q == j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : ScoreSingle(i,j,p,q)),
                                   EncodeTraceback(TB_FC_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q));
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    int bp = -1, bq = -1;
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCv[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            RealT score = (p == i && q == j) ?
                                (score_helix + FCptr[q]) :
                                (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                            if (score > best_v)
                            {
                                best_v = score;
                                bp = p;
                                bq = q;
                            }
                        }
                    }
                    
                    if (bp != -1 && bq != -1)
                        best_t = EncodeTraceback(TB_FC_SINGLE,(bp-i)*(C_MAX_SINGLE_LENGTH+1)+j-bq);
                }
#endif
                
                // compute MAX (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                UPDATE_MAX(best_v, best_t,
                           FM2v + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase(), 
                           EncodeTraceback(TB_FC_BIFURCATION,FM2t));
                
                FCv[offset[i]+j] = best_v;
                FCt[offset[i]+j] = best_t;
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = MAX [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                {
                    UPDATE_MAX(best_v, best_t, 
                               FCv[offset[i+1]+j-1] + ScoreJunctionA(j,i) +
                               ScoreMultiPaired() + ScoreBasePair(i+1,j), 
                               EncodeTraceback(TB_FM1_PAIRED,0));
                }
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                {
                    UPDATE_MAX(best_v, best_t,
                               FM1v[offset[i+1]+j] + ScoreMultiUnpaired(i+1),
                               EncodeTraceback(TB_FM1_UNPAIRED,0));
                }
                
                FM1v[offset[i]+j] = best_v;
                FM1t[offset[i]+j] = best_t;
            }
            
#if CANDIDATE_LIST
            
            // If there exists some i <= k < j for which
            //   FM1[i,k] + FM[k,j] >= FM1[i,j]
            // then for all j' > j, we know that
            //   FM1[i,k] + FM[k,j'] >= FM1[i,j] + FM[j,j'].
            // since 
            //   FM[k,j'] >= FM[k,j] + FM[j,j'].
            //
            // From this, it follows that we only need to consider
            // j as a candidate partition point for future j' values
            // only if FM1[i,j] > FM1[i,k] + FM[k,j] for all k.
            
            if (FM1v[offset[i]+j] > FM2v)
                candidates.push_back(j);
#endif
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = MAX [MAX (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute MAX (i<k<j : FM1[i,k] + FM[k,j])
                
                UPDATE_MAX(best_v, best_t, FM2v, EncodeTraceback(TB_FM_BIFURCATION,FM2t));
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                {
                    UPDATE_MAX(best_v, best_t,
                               FMv[offset[i]+j-1] + ScoreMultiUnpaired(j), 
                               EncodeTraceback(TB_FM_UNPAIRED,0));
                }
                
                // compute FM1[i,j]
                
                UPDATE_MAX(best_v, best_t, FM1v[offset[i]+j], EncodeTraceback(TB_FM_FM1,0));
                
                FMv[offset[i]+j] = best_v;
                FMt[offset[i]+j] = best_t;
            }
        }
    }
    
    F5v[0] = RealT(0);
    F5t[0] = EncodeTraceback(TB_F5_ZERO,0);
    for (int j = 1; j <= L; j++)
    {
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = MAX [F5[j-1] + ScoreExternalUnpaired(),
        //              MAX (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT best_v = RealT(NEG_INF);
        int best_t = -1;
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
        {
            UPDATE_MAX(best_v, best_t, 
                       F5v[j-1] + ScoreExternalUnpaired(j),
                       EncodeTraceback(TB_F5_UNPAIRED,0));
        }
        
        // compute MAX (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j])
            {
                UPDATE_MAX(best_v, best_t,
                           F5v[k] + FCv[offset[k+1]+j-1] + ScoreExternalPaired() +
                           ScoreBasePair(k+1,j) + ScoreJunctionA(j,k),
                           EncodeTraceback(TB_F5_BIFURCATION,k));
            }
        }
        
        F5v[j] = best_v;
        F5t[j] = best_t;
    }

#if SHOW_TIMINGS
    std::cerr << "Viterbi score: " << F5v[L] << " (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif
    
#if CANDIDATE_LIST
    //std::cerr << "Candidates: " << candidates_seen << "/" << candidates_possible << " = " << double(candidates_seen)/candidates_possible << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::GetViterbiScore()
//
// Return Viterbi score for a sequence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::GetViterbiScore() const
{
    return F5v[L];
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::PredictPairingsViterbi()
// 
// Use Viterbi decoding to predict pairings.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<int> InferenceEngine<RealT>::PredictPairingsViterbi() const
{
    std::vector<int> solution(L+1,SStruct::UNPAIRED);
    solution[0] = SStruct::UNKNOWN;
    //return solution;
    
    std::queue<triple<const int *,int,int> > traceback_queue;
    traceback_queue.push(make_triple(&F5t[0], 0, L));
    
    while (!traceback_queue.empty())
    {
        triple<const int *,int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int *V = t.first;
        const int i = t.second;
        const int j = t.third;
        
        std::pair<int,int> traceback = DecodeTraceback(V == &F5t[0] ? V[j] : V[offset[i]+j]);
        
        //std::cerr << (V == FCt ? "FC " : V == FMt ? "FM " : V == FM1t ? "FM1 " : "F5 ");
        //std::cerr << i << " " << j << ": " << traceback.first << " " << traceback.second << std::endl;
        
        switch (traceback.first)
        {
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            case TB_FN_HAIRPIN: 
                break;
            case TB_FN_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                solution[p+1] = q;
                solution[q] = p+1;
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FN_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FE_STACKING: 
            {
                solution[i+1] = j;
                solution[j] = i+1;
                traceback_queue.push(make_triple(&FEt[0], i+1, j-1));
            }
            break;
            case TB_FE_FN: 
            {
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_FN:
            {
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_HELIX:
            {
                const int m = traceback.second;
                for (int k = 2; k <= m; k++)
                {
                    solution[i+k-1] = j-k+2;
                    solution[j-k+2] = i+k-1;
                }
                traceback_queue.push(make_triple(&FNt[0], i+m-1, j-m+1));
            }
            break;
            case TB_FC_FE:
            {
                const int m = D_MAX_HELIX_LENGTH;
                for (int k = 2; k <= m; k++)
                {
                    solution[i+k-1] = j-k+2;
                    solution[j-k+2] = i+k-1;
                }
                traceback_queue.push(make_triple(&FEt[0], i+m-1, j-m+1));
            }
            break;
#else
            case TB_FC_HAIRPIN: 
                break;
            case TB_FC_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                solution[p+1] = q;
                solution[q] = p+1;
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FC_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
#endif
            case TB_FM1_PAIRED:
            {
                solution[i+1] = j;
                solution[j] = i+1;
                traceback_queue.push(make_triple(&FCt[0], i+1, j-1));
            }
      break;
            case TB_FM1_UNPAIRED:
            {
                traceback_queue.push(make_triple(&FM1t[0], i+1, j));
            }
            break;
            case TB_FM_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FM_UNPAIRED:
            {
                traceback_queue.push(make_triple(&FMt[0], i, j-1));
            }
            break;
            case TB_FM_FM1: 
            {
                traceback_queue.push(make_triple(&FM1t[0], i, j));
            }
            break;
            case TB_F5_ZERO:
                break;
            case TB_F5_UNPAIRED:
            {
                traceback_queue.push(make_triple(&F5t[0], 0, j-1));
            }
            break;
            case TB_F5_BIFURCATION:
            {
                const int k = traceback.second;
                solution[k+1] = j;
                solution[j] = k+1;
                traceback_queue.push(make_triple(&F5t[0], 0, k));
                traceback_queue.push(make_triple(&FCt[0], k+1, j-1));
            }
            break;
            default:
                Assert(false, "Bad traceback.");
        }
    }
    
    return solution;
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeViterbiFeatureCounts()
// 
// Use feature counts from Viterbi decoding.
//////////////////////////////////////////////////////////////////////

template<class RealT>
typename InferenceEngine<RealT>::CntPtr
InferenceEngine<RealT>::ComputeViterbiFeatureCounts()
{
    std::queue<triple<int *,int,int> > traceback_queue;
    traceback_queue.push(make_triple(&F5t[0], 0, L));

    ClearCounts();
    parameter_count.reset(new ParameterHash<uint>());
    
    while (!traceback_queue.empty())
    {
        triple<int *,int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int *V = t.first;
        const int i = t.second;
        const int j = t.third;
        
        std::pair<int,int> traceback = DecodeTraceback (V == &F5t[0] ? V[j] : V[offset[i]+j]);
        
        switch (traceback.first)
        {
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            case TB_FN_HAIRPIN: 
                CountHairpin(i,j,1);
                break;
            case TB_FN_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                CountSingle(i,j,p,q,1);
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FN_BIFURCATION:
            {
                const int k = traceback.second;
                CountJunctionA(i,j,1);
                CountMultiPaired(1);
                CountMultiBase(1);
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FE_STACKING: 
            {
                CountBasePair(i+1,j,1);
                CountHelixStacking(i,j+1,1);
                traceback_queue.push(make_triple(&FEt[0], i+1, j-1));
            }
            break;
            case TB_FE_FN: 
            {
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_FN:
            {
                CountIsolated(1);
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_HELIX:
            {
                const int m = traceback.second;
                CountHelix(i-1,j+1,m,1);
                traceback_queue.push(make_triple(&FNt[0], i+m-1, j-m+1));
            }
            break;
            case TB_FC_FE:
            {
                const int m = D_MAX_HELIX_LENGTH;
                CountHelix(i-1,j+1,m,1);
                traceback_queue.push(make_triple(&FEt[0], i+m-1, j-m+1));
            }
            break;
#else
            case TB_FC_HAIRPIN: 
                CountHairpin(i,j,1);
                break;
            case TB_FC_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);

                if (p == i && q == j)
                {
                    CountBasePair(i+1,j,1);
                    CountHelixStacking(i,j+1,1);
                }
                else
                {
                    CountSingle(i,j,p,q,1);
                }
                
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FC_BIFURCATION:
            {
                const int k = traceback.second;
                CountJunctionA(i,j,1);
                CountMultiPaired(1);
                CountMultiBase(1);
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
#endif
            case TB_FM1_PAIRED:
            {
                CountJunctionA(j,i,1);
                CountMultiPaired(1);
                CountBasePair(i+1,j,1);
                traceback_queue.push(make_triple(&FCt[0], i+1, j-1));
            }
            break;
            case TB_FM1_UNPAIRED:
            {
                CountMultiUnpaired(i+1,1);
                traceback_queue.push(make_triple(&FM1t[0], i+1, j));
            }
            break;
            case TB_FM_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FM_UNPAIRED:
            {
                CountMultiUnpaired(j,1);
                traceback_queue.push(make_triple(&FMt[0], i, j-1));
            }
            break;
            case TB_FM_FM1: 
                traceback_queue.push(make_triple(&FM1t[0], i, j));
                break;
            case TB_F5_ZERO:
                break;
            case TB_F5_UNPAIRED:
                CountExternalUnpaired(j,1);
                traceback_queue.push(make_triple(&F5t[0], 0, j-1));
                break;
            case TB_F5_BIFURCATION:
            {
                const int k = traceback.second;
                CountExternalPaired(1);
                CountBasePair(k+1,j,1);
                CountJunctionA(j,k,1);
                traceback_queue.push(make_triple(&F5t[0], 0, k));
                traceback_queue.push(make_triple(&FCt[0], k+1, j-1));
            }
            break;
            default:
                Assert(false, "Bad traceback.");
        }
    }

    FinalizeCounts();
    return std::move(parameter_count);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeInside()
//
// Run inside algorithm.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputeInside()
{
    InitializeCache();
        
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif
    
    // initialization

    F5i.clear(); F5i.resize(L+1, RealT(NEG_INF));
    FCi.clear(); FCi.resize(SIZE, RealT(NEG_INF));
    FMi.clear(); FMi.resize(SIZE, RealT(NEG_INF));
    FM1i.clear(); FM1i.resize(SIZE, RealT(NEG_INF));
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    FEi.clear(); FEi.resize(SIZE, RealT(NEG_INF));
    FNi.clear(); FNi.resize(SIZE, RealT(NEG_INF));
#endif

    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i = RealT(NEG_INF);
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i, FM1i[offset[i]+k] + FMi[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i[offset[i]+i+1]);
                const RealT *p2 = &(FMi[offset[i+1]+j]);
                for (int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
            
#endif
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR

            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    Fast_LogPlusEquals(sum_i, ScoreHairpin(i,j));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        Fast_LogPlusEquals(sum_i, ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]);
                    }
                }
                
#else

                if (i+2 <= j)                
                {
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT score = (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) +
                                           ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                            Fast_LogPlusEquals(sum_i, score);
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(sum_i, FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
                FNi[offset[i]+j] = sum_i;
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    Fast_LogPlusEquals(sum_i, ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FEi[offset[i+1]+j-1]);
                }
                
                // compute FN(i,j)

                Fast_LogPlusEquals(sum_i, FNi[offset[i]+j]);
                
                FEi[offset[i]+j] = sum_i;
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreIsolated() + FN(i,j)
                
                Fast_LogPlusEquals(sum_i, ScoreIsolated() + FNi[offset[i]+j]);
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    Fast_LogPlusEquals(sum_i, ScoreHelix(i-1,j+1,k) + FNi[offset[i+k-1]+j-k+1]);
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        Fast_LogPlusEquals(sum_i, ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]);
                }
                
                FCi[offset[i]+j] = sum_i;
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    Fast_LogPlusEquals(sum_i, ScoreHairpin(i,j));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        Fast_LogPlusEquals(sum_i,
                                           FCi[offset[p+1]+q-1] +
                                           (p == i && q == j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : ScoreSingle(i,j,p,q)));
                    }
                }
                
#else

                {
                    RealT score_helix = (i+2 <= j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            RealT score = (p == i && q == j) ?
                                (score_helix + FCptr[q]) :
                                (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) +
                                 ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                            Fast_LogPlusEquals(sum_i, score);
                        }
                    }
                }
#endif

                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(sum_i, FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
                FCi[offset[i]+j] = sum_i;
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    Fast_LogPlusEquals(sum_i, FCi[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j));
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                    Fast_LogPlusEquals(sum_i, FM1i[offset[i+1]+j] + ScoreMultiUnpaired(i+1));
                
                FM1i[offset[i]+j] = sum_i;
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j])
                
                Fast_LogPlusEquals(sum_i, FM2i);
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                    Fast_LogPlusEquals(sum_i, FMi[offset[i]+j-1] + ScoreMultiUnpaired(j));
                
                // compute FM1[i,j]
                
                Fast_LogPlusEquals(sum_i, FM1i[offset[i]+j]);
                
                FMi[offset[i]+j] = sum_i;
            }
        }
    }
    
    F5i[0] = RealT(0);
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT sum_i = RealT(NEG_INF);
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            Fast_LogPlusEquals(sum_i, F5i[j-1] + ScoreExternalUnpaired(j));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
            if (allow_paired[offset[k+1]+j])
                Fast_LogPlusEquals(sum_i, F5i[k] + FCi[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k));
        
        F5i[j] = sum_i;
    }

#if SHOW_TIMINGS
    std::cerr << "Inside score: " << F5i[L] << " (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeOutside()
//
// Run outside algorithm.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputeOutside()
{
    InitializeCache();
    
#if SHOW_TIMINGS    
    double starting_time = GetSystemTime();
#endif
    
    // initialization
    
    F5o.clear(); F5o.resize(L+1, RealT(NEG_INF));
    FCo.clear(); FCo.resize(SIZE, RealT(NEG_INF));
    FMo.clear(); FMo.resize(SIZE, RealT(NEG_INF));
    FM1o.clear(); FM1o.resize(SIZE, RealT(NEG_INF));
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    FEo.clear(); FEo.resize(SIZE, RealT(NEG_INF));
    FNo.clear(); FNo.resize(SIZE, RealT(NEG_INF));
#endif
    
    F5o[L] = RealT(0);  
    for (int j = L; j >= 1; j--)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            Fast_LogPlusEquals(F5o[j-1], F5o[j] + ScoreExternalUnpaired(j));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        {
            for (int k = 0; k < j; k++)
            {
                if (allow_paired[offset[k+1]+j])
                {
                    RealT temp = F5o[j] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k);
                    Fast_LogPlusEquals(F5o[k], temp + FCi[offset[k+1]+j-1]);
                    Fast_LogPlusEquals(FCo[offset[k+1]+j-1], temp + F5i[k]);
                }
            }
        }
    }
    
    for (int i = 0; i <= L; i++)
    {
        for (int j = L; j >= i; j--)
        {
            RealT FM2o = RealT(NEG_INF);
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j])
                
                Fast_LogPlusEquals(FM2o, FMo[offset[i]+j]);
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                    Fast_LogPlusEquals(FMo[offset[i]+j-1], FMo[offset[i]+j] + ScoreMultiUnpaired(j));
                
                // compute FM1[i,j]
                
                Fast_LogPlusEquals(FM1o[offset[i]+j], FMo[offset[i]+j]);
            }
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    Fast_LogPlusEquals(FCo[offset[i+1]+j-1], FM1o[offset[i]+j] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j));
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                    Fast_LogPlusEquals(FM1o[offset[i+1]+j], FM1o[offset[i]+j] + ScoreMultiUnpaired(i+1));
                
            }
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreIsolated() + FN(i,j)
                
                Fast_LogPlusEquals(FNo[offset[i]+j], ScoreIsolated() + FCo[offset[i]+j]);
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    Fast_LogPlusEquals(FNo[offset[i+k-1]+j-k+1], ScoreHelix(i-1,j+1,k) + FCo[offset[i]+j]);
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        Fast_LogPlusEquals(FEo[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1],
                                           ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) + FCo[offset[i]+j]);
                }
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    Fast_LogPlusEquals(FEo[offset[i+1]+j-1], FEo[offset[i]+j] + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1));
                }
                
                // compute FN(i,j)
                
                Fast_LogPlusEquals(FNo[offset[i]+j], FEo[offset[i]+j]);
            }
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreHairpin(i,j) -- do nothing
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                {
                    RealT temp = FNo[offset[i]+j];
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            Fast_LogPlusEquals(FCo[offset[p+1]+q-1], temp + ScoreSingle(i,j,p,q));
                        }
                    }
                }
#else
                
                {
                    RealT score_other = FNo[offset[i]+j] + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        RealT *FCptr = &(FCo[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            Fast_LogPlusEquals(FCptr[q], score_other + cache_score_single[p-i][j-q].first + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(FM2o, FNo[offset[i]+j] + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                // compute ScoreHairpin(i,j) -- do nothing
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])

#if !FAST_SINGLE_BRANCH_LOOPS
                {
                    RealT temp = FCo[offset[i]+j];
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            Fast_LogPlusEquals(FCo[offset[p+1]+q-1],
                                               temp + (p == i && q == j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : ScoreSingle(i,j,p,q)));
                        }
                    }
                }
#else
                
                {
                    RealT score_helix = (i+2 <= j ? FCo[offset[i]+j] + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = FCo[offset[i]+j] + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        RealT *FCptr = &(FCo[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            Fast_LogPlusEquals(FCptr[q], 
                                               (p == i && q == j) ?
                                               score_helix :
                                               score_other + cache_score_single[p-i][j-q].first + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(FM2o, FCo[offset[i]+j] + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
            }
            
#endif
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
            {
                Fast_LogPlusEquals(FM1o[offset[i]+k], FM2o + FMi[offset[k]+j]);
                Fast_LogPlusEquals(FMo[offset[k]+j], FM2o + FM1i[offset[i]+k]);
            }

#else
            if (i+2 <= j)
            {
                RealT *p1i = &(FM1i[offset[i]+i+1]);
                RealT *p2i = &(FMi[offset[i+1]+j]);
                RealT *p1o = &(FM1o[offset[i]+i+1]);
                RealT *p2o = &(FMo[offset[i+1]+j]);
                for (int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(*p1o, FM2o + *p2i);
                    Fast_LogPlusEquals(*p2o, FM2o + *p1i);
                    ++p1i;
                    ++p1o;
                    p2i += L-k;
                    p2o += L-k;
                }
            }
            
#endif
        }
    }
    
#if SHOW_TIMINGS
    std::cerr << "Outside score: " << F5o[0] << " (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeLogPartitionCoefficient()
//
// Return partition coefficient.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ComputeLogPartitionCoefficient() const
{
    // NOTE: This should be equal to F5o[0]. 
    
    return F5i[L];
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeFeatureCountExpectations()
// 
// Combine the results of the inside and outside algorithms
// in order to compute feature count expectations.
//////////////////////////////////////////////////////////////////////

template<class RealT>
typename InferenceEngine<RealT>::CntPtr 
InferenceEngine<RealT>::ComputeFeatureCountExpectations()
{
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif

    //std::cerr << "Inside score: " << F5i[L].GetLogRepresentation() << std::endl;
    //std::cerr << "Outside score: " << F5o[0].GetLogRepresentation() << std::endl;
    
    const RealT Z = ComputeLogPartitionCoefficient();

    ClearCounts();
    parameter_count.reset(new ParameterHash<uint>());
    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {

            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i = RealT(NEG_INF);
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i, FM1i[offset[i]+k] + FMi[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i[offset[i]+i+1]);
                const RealT *p2 = &(FMi[offset[i+1]+j]);
                for (int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
            
#endif
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FNo[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpin(i,j,Fast_Exp(outside + ScoreHairpin(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        CountSingle(i,j,p,q,Fast_Exp(outside + ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]));
                    }
                }
                
#else
                
                {
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT value = Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            cache_score_single[p-i][j-q].second += value;
                            CountBasePair(p+1,q,value);
                            CountJunctionB(i,j,value);
                            CountJunctionB(q,p,value);
                            CountSingleNucleotides(i,j,p,q,value);
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                {
                    RealT value = Fast_Exp(outside + FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                    CountJunctionA(i,j,value);
                    CountMultiPaired(value);
                    CountMultiBase(value);
                }
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FEo[offset[i]+j] - Z;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    RealT value = Fast_Exp(outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FEi[offset[i+1]+j-1]);
                    CountBasePair(i+1,j,value);
                    CountHelixStacking(i,j+1,value);
                }
                
                // compute FN(i,j) -- do nothing
                
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FCo[offset[i]+j] - Z;
                
                // compute ScoreIsolated() + FN(i,j)
                
                CountIsolated(Fast_Exp(outside + ScoreIsolated() + FNi[offset[i]+j]));
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    CountHelix(i-1,j+1,k,Fast_Exp(outside + ScoreHelix(i-1,j+1,k) + FNi[offset[i+k-1]+j-k+1]));
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        CountHelix(i-1,j+1,D_MAX_HELIX_LENGTH,
                                   Fast_Exp(outside + ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]));
                }
            }

#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FCo[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpin(i,j,Fast_Exp(outside + ScoreHairpin(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        if (p == i && q == j)
                        {
                            RealT value = Fast_Exp(outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FCi[offset[p+1]+q-1]);
                            CountBasePair(i+1,j,value);
                            CountHelixStacking(i,j+1,value);
                        }
                        else
                        {
                            CountSingle(i,j,p,q,Fast_Exp(outside + ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]));
                        }
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            if (p == i && q == j)
                            {
                                RealT value = Fast_Exp(score_helix + FCptr[q]);
                                cache_score_single[0][0].second += value;
                                CountBasePair(i+1,j,value);
                                CountHelixStacking(i,j+1,value);
                            }
                            else
                            {
                                RealT value = Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + 
                                                       ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                                cache_score_single[p-i][j-q].second += value;
                                CountBasePair(p+1,q,value);
                                CountJunctionB(i,j,value);
                                CountJunctionB(q,p,value);
                                CountSingleNucleotides(i,j,p,q,value);
                            }
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                {
                    RealT value = Fast_Exp(outside + FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                    CountJunctionA(i,j,value);
                    CountMultiPaired(value);
                    CountMultiBase(value);
                }
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                {
                    RealT value = Fast_Exp(FM1o[offset[i]+j] + FCi[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j) - Z);
                    CountJunctionA(j,i,value);
                    CountMultiPaired(value);
                    CountBasePair(i+1,j,value);
                }
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                {
                    CountMultiUnpaired(i+1,Fast_Exp(FM1o[offset[i]+j] + FM1i[offset[i+1]+j] + ScoreMultiUnpaired(i+1) - Z));
                }
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j]) -- do nothing
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                    CountMultiUnpaired(j,Fast_Exp(FMo[offset[i]+j] + FMi[offset[i]+j-1] + ScoreMultiUnpaired(j) - Z));
                
                // compute FM1[i,j] -- do nothing
            }
        }
    }
    
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT outside = F5o[j] - Z;
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            CountExternalUnpaired(j,Fast_Exp(outside + F5i[j-1] + ScoreExternalUnpaired(j)));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j])
            {
                RealT value = Fast_Exp(outside + F5i[k] + FCi[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k));
                CountExternalPaired(value);
                CountBasePair(k+1,j,value);
                CountJunctionA(j,k,value);
            }      
        }
    }
    
    FinalizeCounts();

#if SHOW_TIMINGS
    std::cerr << "Feature expectations (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif

    return std::move(parameter_count);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputePosterior()
// 
// Combine the results of the inside and outside algorithms
// in order to compute posterior probabilities of base pairing.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputePosterior()
{ 
    posterior.clear();
    posterior.resize(SIZE, RealT(0));
    
    //double starting_time = GetSystemTime();

    const RealT Z = ComputeLogPartitionCoefficient();
    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i = RealT(NEG_INF);
            
#if SIMPLE_FM2
      
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i, FM1i[offset[i]+k] + FMi[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i[offset[i]+i+1]);
                const RealT *p2 = &(FMi[offset[i+1]+j]);
                for (int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
      
#endif

#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FNo[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpin(i,j,Fast_Exp(outside + ScoreHairpin(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]);
                    }
                }
                
#else
                
                {
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            posterior[offset[p+1]+q] += 
                                Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                        }
                    }
                }
#endif

                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c) -- do nothing
                
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FEo[offset[i]+j] - Z;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                    posterior[offset[i]+j] += Fast_Exp(outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FEi[offset[i+1]+j-1]);
                
                // compute FN(i,j) -- do nothing
                
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FCo[offset[i]+j] - Z;
                
                // compute ScoreIsolated() + FN(i,j) -- do nothing
                
                CountIsolated(Fast_Exp(outside + ScoreIsolated() + FNi[offset[i]+j]));
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    RealT value = Fast_Exp(outside + ScoreHelix(i-1,j+1,k) + FNi[offset[i+k-1]+j-k+1]);
                    for (int p = 1; p < k; p++)
                        posterior[offset[i+p]+j-p+1] += value;
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2]) {
                        RealT value = Fast_Exp(outside + ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]);
                        
                        for (int k = 1; k < D_MAX_HELIX_LENGTH; k++)
                            posterior[offset[i+k]+j-k+1] += value;
                    }
                }
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FCo[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpin(i,j,Fast_Exp(outside + ScoreHairpin(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        if (p == i && q == j)
                        {
                            posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FCi[offset[p+1]+q-1]);
                        }
                        else
                        {
                            posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]);
                        }
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            posterior[offset[p+1]+q] +=
                                Fast_Exp(p == i && q == j ?
                                         score_helix + FCptr[q] :
                                         score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + 
                                         ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                        }
                    }
                }
                
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c) -- do nothing
                
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // Compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    posterior[offset[i+1]+j] += Fast_Exp(FM1o[offset[i]+j] + FCi[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j) - Z);
                
                // Compute FM1[i+1,j] + b -- do nothing
                
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            // Compute SUM (i<k<j : FM1[i,k] + FM[k,j]) -- do nothing
            
            // Compute FM[i,j-1] + b -- do nothing
            
            // Compute FM1[i,j] -- do nothing
        }
    }
    
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT outside = F5o[j] - Z;
        
        // compute F5[j-1] + ScoreExternalUnpaired() -- do nothing
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j])
                posterior[offset[k+1]+j] += Fast_Exp(outside + F5i[k] + FCi[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k));
        }
    }

    for (int i = 1; i <= L; i++)
    {
        for (int j = i+1; j <= L; j++)
        {
            posterior[offset[i]+j] = Clip(posterior[offset[i]+j], RealT(0), RealT(1));
        }
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::PredictPairingsPosterior()
//
// Use posterior decoding to predict pairings.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<int> InferenceEngine<RealT>::PredictPairingsPosterior(const RealT gamma) const
{
    Assert(gamma > 0, "Non-negative gamma expected.");
    
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif
    RealT* unpaired_posterior  = new RealT[L+1];
    RealT* score               = new RealT[SIZE];
    int* traceback             = new int[SIZE];
    
    // compute the scores for unpaired nucleotides
    
    for (int i = 1; i <= L; i++)
    {
        unpaired_posterior[i] = RealT(1);
        for (int j = 1; j < i; j++) unpaired_posterior[i] -= posterior[offset[j]+i];
        for (int j = i+1; j <= L; j++) unpaired_posterior[i] -= posterior[offset[i]+j];
    }
    
    for (int i = 1; i <= L; i++) unpaired_posterior[i] /= 2 * gamma;
    
    // initialize matrices
    
    std::fill(score, score+SIZE, RealT(-1.0));
    std::fill(traceback, traceback+SIZE, -1);
    
    // dynamic programming
    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            RealT &this_score = score[offset[i]+j];
            int &this_traceback = traceback[offset[i]+j];
            
            if (i == j)
            {
                UPDATE_MAX(this_score, this_traceback, RealT(0), 0);
            }
            else
            {
                if (allow_unpaired_position[i+1])
                    UPDATE_MAX(this_score, this_traceback, unpaired_posterior[i+1] + score[offset[i+1]+j], 1);
                if (allow_unpaired_position[j])
                    UPDATE_MAX(this_score, this_traceback, unpaired_posterior[j] + score[offset[i]+j-1], 2);
                if (i+2 <= j)
                { 
                    if (allow_paired[offset[i+1]+j])
                        UPDATE_MAX(this_score, this_traceback, posterior[offset[i+1]+j] + score[offset[i+1]+j-1], 3);
                    
#if SIMPLE_FM2
                    
                    for (int k = i+1; k < j; k++)
                        UPDATE_MAX(this_score, this_traceback, score[offset[i]+k] + score[offset[k]+j], k+4);	
                    
#else
                    
                    RealT *p1 = &(score[offset[i]+i+1]);
                    RealT *p2 = &(score[offset[i+1]+j]);
                    for (int k = i+1; k < j; k++)
                    {
                        UPDATE_MAX(this_score, this_traceback, (*p1) + (*p2), k+4);
                        ++p1;
                        p2 += L-k;
                    }
                    
#endif
                }
            }
        }
    }

#if SHOW_TIMINGS
    std::cerr << "Time: " << GetSystemTime() - starting_time << std::endl;
#endif
    
    // perform traceback
    
    std::vector<int> solution(L+1,SStruct::UNPAIRED);
    solution[0] = SStruct::UNKNOWN;
    
    std::queue<std::pair<int,int> > traceback_queue;
    traceback_queue.push(std::make_pair(0, L));
    
    while (!traceback_queue.empty())
    {
        std::pair<int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int i = t.first;
        const int j = t.second;
        
        switch (traceback[offset[i]+j])
        {
            case -1:
                Assert(false, "Should not get here.");
                break;
            case 0: 
                break;
            case 1: 
                traceback_queue.push(std::make_pair(i+1,j));
                break;
            case 2: 
                traceback_queue.push(std::make_pair(i,j-1));
                break;
            case 3:
                solution[i+1] = j;
                solution[j] = i+1;
                traceback_queue.push(std::make_pair(i+1,j-1));
                break;
            default:
            {
                const int k = traceback[offset[i]+j] - 4;
                traceback_queue.push(std::make_pair(i,k));
                traceback_queue.push(std::make_pair(k,j));
            }
            break;       
        }
    }
    
    return solution;
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::GetPosterior()
//
// Return posterior probability matrix, thresholded.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT *InferenceEngine<RealT>::GetPosterior(const RealT posterior_cutoff) const
{
    RealT *ret = new RealT[SIZE];
    for (int i = 0; i < SIZE; i++)
        ret[i] = (posterior[i] >= posterior_cutoff ? posterior[i] : RealT(0));
    return ret;
}

template 
class InferenceEngine<double>;

// Local Variables:
// mode: C++
// c-basic-offset: 4
// End:
