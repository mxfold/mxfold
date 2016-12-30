#ifndef PARAMETERHASH_HPP
#define PARAMETERHASH_HPP

#include "Config.hpp"
#include "Utilities.hpp"
#include <string>
#include <unordered_map>
#include <array>

template < class ValueT >
class ParameterHash
{
public:
  //typedef typename ValueT ValueT;

public:
  ParameterHash();
  ~ParameterHash() {}

  void initialize_char_mapping(const std::string& alphabet);
  size_t alphabet_size() const { return alphabet_.size(); }
  BYTE char_mapping(BYTE c) const { return char_mapping_[c]; }

  void ReadFromFile(const std::string& filename);
  void WriteToFile(const std::string& filename) const;

  // access to parameters
#if PARAMS_BASE_PAIR
  ValueT  base_pair(BYTE i, BYTE j) const;
  ValueT& base_pair(BYTE i, BYTE j);
#endif
#if PARAMS_BASE_PAIR_DIST
  ValueT  base_pair_dist_at_least(uint l) const;
  ValueT& base_pair_dist_at_least(uint l);
#endif
#if PARAMS_TERMINAL_MISMATCH
  ValueT  terminal_mismatch(BYTE i1, BYTE j1, BYTE i2, BYTE j2) const;
  ValueT& terminal_mismatch(BYTE i1, BYTE j1, BYTE i2, BYTE j2);  
#endif
#if PARAMS_HAIRPIN_LENGTH
  ValueT  hairpin_length_at_least(uint l) const;
  ValueT& hairpin_length_at_least(uint l);
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
  ValueT  hairpin_3_nucleotides(BYTE i1, BYTE i2, BYTE i3) const;
  ValueT& hairpin_3_nucleotides(BYTE i1, BYTE i2, BYTE i3);  
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
  ValueT  hairpin_4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4) const;
  ValueT& hairpin_4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4);  
#endif
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
  ValueT  hairpin_5_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5) const;
  ValueT& hairpin_5_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5);  
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
  ValueT  hairpin_6_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6) const;
  ValueT& hairpin_6_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6);  
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
  ValueT  hairpin_7_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7) const;
  ValueT& hairpin_7_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7);  
#endif
#if PARAMS_HELIX_LENGTH
  ValueT  helix_length_at_least(uint l) const;
  ValueT& helix_length_at_least(uint l);
#endif
#if PARAMS_ISOLATED_BASE_PAIR
  ValueT  isolated_base_pair() const;
  ValueT& isolated_base_pair();
#endif
#if PARAMS_INTERNAL_EXPLICIT
  ValueT  internal_explicit(uint i, uint j) const;
  ValueT& internal_explicit(uint i, uint j);  
#endif
#if PARAMS_BULGE_LENGTH
  ValueT  bulge_length_at_least(uint l) const;
  ValueT& bulge_length_at_least(uint l);
#endif
#if PARAMS_INTERNAL_LENGTH
  ValueT  internal_length_at_least(uint l) const;
  ValueT& internal_length_at_least(uint l);
#endif
#if PARAMS_INTERNAL_SYMMETRY
  ValueT  internal_symmetric_length_at_least(uint l) const;
  ValueT& internal_symmetric_length_at_least(uint l);
#endif
#if PARAMS_INTERNAL_ASYMMETRY
  ValueT  internal_asymmetry_at_least(uint l) const;
  ValueT& internal_asymmetry_at_least(uint l);
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
  ValueT  bulge_0x1_nucleotides(BYTE i1) const;
  ValueT& bulge_0x1_nucleotides(BYTE i1);
  ValueT  bulge_1x0_nucleotides(BYTE i1) const;
  ValueT& bulge_1x0_nucleotides(BYTE i1);
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
  ValueT  bulge_0x2_nucleotides(BYTE i1, BYTE i2) const;
  ValueT& bulge_0x2_nucleotides(BYTE i1, BYTE i2);
  ValueT  bulge_2x0_nucleotides(BYTE i1, BYTE i2) const;
  ValueT& bulge_2x0_nucleotides(BYTE i1, BYTE i2);
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
  ValueT  bulge_0x3_nucleotides(BYTE i1, BYTE i2, BYTE i3) const;
  ValueT& bulge_0x3_nucleotides(BYTE i1, BYTE i2, BYTE i3);
  ValueT  bulge_3x0_nucleotides(BYTE i1, BYTE i2, BYTE i3) const;
  ValueT& bulge_3x0_nucleotides(BYTE i1, BYTE i2, BYTE i3);
#endif
#if PARAMS_BULGE_0x4_NUCLEOTIDES
  ValueT  bulge_0x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4) const;
  ValueT& bulge_0x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4);
  ValueT  bulge_4x0_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4) const;
  ValueT& bulge_4x0_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4);
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
  ValueT  bulge_0x5_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5) const;
  ValueT& bulge_0x5_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5);
  ValueT  bulge_5x0_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5) const;
  ValueT& bulge_5x0_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5);
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
  ValueT  bulge_0x6_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6) const;
  ValueT& bulge_0x6_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6);
  ValueT  bulge_6x0_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6) const;
  ValueT& bulge_6x0_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6);
#endif
#if PARAMS_BULGE_0x7_NUCLEOTIDES
  ValueT  bulge_0x7_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7) const;
  ValueT& bulge_0x7_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7);
  ValueT  bulge_7x0_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7) const;
  ValueT& bulge_7x0_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7);
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
  ValueT  internal_1x1_nucleotides(BYTE i1, BYTE i2) const;
  ValueT& internal_1x1_nucleotides(BYTE i1, BYTE i2);
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
  ValueT  internal_1x2_nucleotides(BYTE i1, BYTE i2, BYTE i3) const;
  ValueT& internal_1x2_nucleotides(BYTE i1, BYTE i2, BYTE i3);
  ValueT  internal_2x1_nucleotides(BYTE i1, BYTE i2, BYTE i3) const;
  ValueT& internal_2x1_nucleotides(BYTE i1, BYTE i2, BYTE i3);
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
  ValueT  internal_2x2_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4) const;
  ValueT& internal_2x2_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4);
#endif
#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
  ValueT  internal_1x3_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4) const;
  ValueT& internal_1x3_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4);
  ValueT  internal_3x1_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4) const;
  ValueT& internal_3x1_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4);
#endif
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
  ValueT  internal_2x3_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5) const;
  ValueT& internal_2x3_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5);
  ValueT  internal_3x2_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5) const;
  ValueT& internal_3x2_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5);
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
  ValueT  internal_3x3_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6) const;
  ValueT& internal_3x3_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6);
#endif

#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
  ValueT  internal_1x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5) const;
  ValueT& internal_1x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5);
  ValueT  internal_4x1_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5) const;
  ValueT& internal_4x1_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5);
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
  ValueT  internal_2x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6) const;
  ValueT& internal_2x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6);
  ValueT  internal_4x2_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6) const;
  ValueT& internal_4x2_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6);
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
  ValueT  internal_3x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7) const;
  ValueT& internal_3x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7);
  ValueT  internal_4x3_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7) const;
  ValueT& internal_4x3_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7);
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
  ValueT  internal_4x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7, BYTE i8) const;
  ValueT& internal_4x4_nucleotides(BYTE i1, BYTE i2, BYTE i3, BYTE i4, BYTE i5, BYTE i6, BYTE i7, BYTE i8);
#endif
#if PARAMS_HELIX_STACKING
  ValueT  helix_stacking(BYTE i1, BYTE j1, BYTE i2, BYTE j2) const;
  ValueT& helix_stacking(BYTE i1, BYTE j1, BYTE i2, BYTE j2);
#endif
#if PARAMS_HELIX_CLOSING
  ValueT  helix_closing(BYTE i, BYTE j) const;
  ValueT& helix_closing(BYTE i, BYTE j);
#endif
#if PARAMS_MULTI_LENGTH
  ValueT  multi_base() const;
  ValueT& multi_base();
  ValueT  multi_unpaired() const;
  ValueT& multi_unpaired();
  ValueT  multi_paired() const;
  ValueT& multi_paired();
#endif
#if PARAMS_DANGLE
  ValueT  dangle_left(BYTE i1, BYTE j1, BYTE i2) const;
  ValueT& dangle_left(BYTE i1, BYTE j1, BYTE i2);
  ValueT  dangle_right(BYTE i1, BYTE j1, BYTE j2) const;
  ValueT& dangle_right(BYTE i1, BYTE j1, BYTE j2);
#endif
#if PARAMS_EXTERNAL_LENGTH
  ValueT  external_unpaired() const;
  ValueT& external_unpaired();
  ValueT  external_paired() const;
  ValueT& external_paired();
#endif

private:
  ValueT  get_by_key(const std::string& key) const;
  ValueT& get_by_key(const std::string& key);

private:
  std::unordered_map<std::string, ValueT> param_;
  std::string alphabet_;
  std::array<BYTE, 256> char_mapping_;
};

#endif  // PARAMETERHASH_HPP

// Local Variables:
// mode: C++
// End:

