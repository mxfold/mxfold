//////////////////////////////////////////////////////////////////////
// SStruct.hpp
//
// This is a class for reading and writing of RNA secondary
// structures.  The file formats supported include
// 
//     (1) BPSEQ
//     (2) FASTA
//     (3) plain text (raw)
//////////////////////////////////////////////////////////////////////

#ifndef SSTRUCT_HPP
#define SSTRUCT_HPP

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// class SStruct
//////////////////////////////////////////////////////////////////////
//template<class RealT>
class SStruct
{
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    std::vector<int> mapping;
    std::vector<float> reactivity_unpair;
    std::vector<float> reactivity_pair;
    int type;
 
    // automatic file format detection
    int AnalyzeFormat(const std::string &filename) const;

    // perform standard character conversions for RNA sequence and structures
    std::string FilterSequence(std::string sequence) const;
    std::string FilterParens(std::string sequence) const;

    // convert a pseudoknot-free parenthesized structure to a mapping and back
    std::vector<int> ConvertParensToMapping(const std::string &parens) const;
    std::string ConvertMappingToParens(const std::vector<int> &mapping) const;

    // check that a (possibly pseudoknotted) mapping is valid
    void ValidateMapping(const std::vector<int> &mapping) const;
    
public:

    // integer constants used to identify nucleotides which are either
    // unpaired or whose pairing is not known
    enum { UNPAIRED = 0, UNKNOWN = -1, PAIRED = -2 };

    // specify the type of the BPSEQ format that contains one of the following information
    enum { NO_REACTIVITY, REACTIVITY_UNPAIRED, REACTIVITY_PAIRED, REACTIVITY_BOTH };

    // constructor and destructor
    SStruct();
    SStruct(const std::string &filename, int type = NO_REACTIVITY);
    SStruct(const SStruct &rhs);
    ~SStruct();

    // load sequence and struture from file
    void Load(const std::string &filename, int type = NO_REACTIVITY);

    // load file of a particular file format
    void LoadFASTA(const std::string &filename);
    void LoadRAW(const std::string &filename);
    void LoadBPSEQ(const std::string &filename, int type = NO_REACTIVITY);

    // assignment operator
    const SStruct& operator=(const SStruct &rhs);

    // check for pseudoknots
    bool ContainsPseudoknots() const;

    // remove noncomplementary base-pairs
    void RemoveNoncomplementaryPairings(const int seq = 0);

    // output in various formats
    void WriteBPSEQ(std::ostream &outfile, const int seq = 0) const;
    void WriteParens(std::ostream &outfile) const;

    // compute alignment percent identity
    double ComputePercentIdentity() const;

    // compute position-based sequence weights
    std::vector<double> ComputePositionBasedSequenceWeights() const;

    // set mapping
    void SetMapping(const std::vector<int> &mapping);

    // discretize reactivity
    void DiscretizeReactivity(double threshold_unpaired, double threshold_paired);

    //////////////////////////////////////////////////////////////////////
    // Getters
    //////////////////////////////////////////////////////////////////////
    
    const std::vector<std::string> &GetNames() const { return names; }
    const std::vector<std::string> &GetSequences() const { return sequences; }
    const std::vector<int> &GetMapping() const { return mapping; }
    int GetLength() const { return int(mapping.size())-1; }
    int GetNumSequences() const { return int(sequences.size()); }
    //get reactivity
    const std::vector<float> &GetReactivityUnpair() const { return reactivity_unpair; }
    const std::vector<float> &GetReactivityPair() const { return reactivity_pair; }
    int GetType() const { return this->type; }
};

#endif

// Local Variables:
// mode: C++
// c-basic-offset: 4
// End:
