#ifndef PARAMETERHASH_HPP
#define PARAMETERHASH_HPP

#include "Config.hpp"
#include "Utilities.hpp"
#include <string>
#include <unordered_map>
#include <iterator>
#include <array>

#define USE_CACHE

template < class ValueT >
class ParameterHash
{
public:
  //typedef typename ValueT ValueT;
  typedef typename std::vector<std::pair<std::string,ValueT>>::iterator iterator;
  typedef typename std::vector<std::pair<std::string,ValueT>>::const_iterator const_iterator;

public:
  ParameterHash();
  ~ParameterHash() { }

  void initialize();
  void initialize_cache();

  bool is_complementary(NUCL x, NUCL y) const;
  bool is_base(NUCL x) const;

  void LoadFromHash(const std::unordered_map<std::string, ValueT>& hash);
  void LoadDefaultComplementary();
  void LoadDefaultNonComplementary();
  void ReadFromFile(const std::string& filename);
  void WriteToFile(const std::string& filename) const;

  ValueT  get_by_key(const std::string& key) const;
  ValueT& get_by_key(const std::string& key);
  uint set_key_value(const std::string& key, ValueT val);
  uint set_key_value(const std::string& key);

  bool is_basepair_feature(const std::string& f) const;
  bool is_context_feature(const std::string& f) const;
  bool is_basepair_context_feature(const std::string& f) const;

  iterator begin()
  {
    return keyval_.begin();
  }

  iterator end()
  {
    return keyval_.end();
  }

  const_iterator begin() const
  {
    return keyval_.begin();
  }

  const_iterator end() const
  {
    return keyval_.end();
  }

  // access to parameters
#if PARAMS_BASE_PAIR
  ValueT  base_pair(NUCL i, NUCL j) const;
  ValueT& base_pair(NUCL i, NUCL j);
  void initialize_cache_base_pair();
#endif
#if PARAMS_BASE_PAIR_DIST
  ValueT  base_pair_dist_at_least(uint l) const;
  ValueT& base_pair_dist_at_least(uint l);
  void initialize_cache_base_pair_dist_at_least();
#endif
#if PARAMS_TERMINAL_MISMATCH
  ValueT  terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const;
  ValueT& terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2);
  void initialize_cache_terminal_mismatch();
#endif
#if PARAMS_HAIRPIN_LENGTH
  ValueT  hairpin_length_at_least(uint l) const;
  ValueT& hairpin_length_at_least(uint l);
  void initialize_cache_hairpin_length_at_least();
#endif

  ValueT  hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l) const;
  ValueT& hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l);

#if PARAMS_HELIX_LENGTH
  ValueT  helix_length_at_least(uint l) const;
  ValueT& helix_length_at_least(uint l);
  void initialize_cache_helix_length_at_least();
#endif
#if PARAMS_ISOLATED_BASE_PAIR
  ValueT  isolated_base_pair() const;
  ValueT& isolated_base_pair();
  void initialize_cache_isolated_base_pair();
#endif
#if PARAMS_INTERNAL_EXPLICIT
  ValueT  internal_explicit(uint i, uint j) const;
  ValueT& internal_explicit(uint i, uint j);
  void initialize_cache_internal_explicit();
#endif
#if PARAMS_BULGE_LENGTH
  ValueT  bulge_length_at_least(uint l) const;
  ValueT& bulge_length_at_least(uint l);
  void initialize_cache_bulge_length_at_least();
#endif
#if PARAMS_INTERNAL_LENGTH
  ValueT  internal_length_at_least(uint l) const;
  ValueT& internal_length_at_least(uint l);
  void initialize_cache_internal_length_at_least();
#endif
#if PARAMS_INTERNAL_SYMMETRY
  ValueT  internal_symmetric_length_at_least(uint l) const;
  ValueT& internal_symmetric_length_at_least(uint l);
  void initialize_cache_internal_symmetric_length_at_least();
#endif
#if PARAMS_INTERNAL_ASYMMETRY
  ValueT  internal_asymmetry_at_least(uint l) const;
  ValueT& internal_asymmetry_at_least(uint l);
  void initialize_cache_internal_asymmetry_at_least();
#endif

  ValueT  internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m) const;
  ValueT& internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m);

#if PARAMS_HELIX_STACKING
  ValueT  helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const;
  ValueT& helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2);
  void initialize_cache_helix_stacking();
#endif
#if PARAMS_HELIX_CLOSING
  ValueT  helix_closing(NUCL i, NUCL j) const;
  ValueT& helix_closing(NUCL i, NUCL j);
  void initialize_cache_helix_closing();
#endif
#if PARAMS_MULTI_LENGTH
  ValueT  multi_base() const;
  ValueT& multi_base();
  void initialize_cache_multi_base();

  ValueT  multi_unpaired() const;
  ValueT& multi_unpaired();
  void initialize_cache_multi_unpaired();

  ValueT  multi_paired() const;
  ValueT& multi_paired();
  void initialize_cache_multi_paired();
#endif
#if PARAMS_DANGLE
  ValueT  dangle_left(NUCL i1, NUCL j1, NUCL i2) const;
  ValueT& dangle_left(NUCL i1, NUCL j1, NUCL i2);
  void initialize_cache_dangle_left();

  ValueT  dangle_right(NUCL i1, NUCL j1, NUCL j2) const;
  ValueT& dangle_right(NUCL i1, NUCL j1, NUCL j2);
  void initialize_cache_dangle_right();
#endif
#if PARAMS_EXTERNAL_LENGTH
  ValueT  external_unpaired() const;
  ValueT& external_unpaired();
  void initialize_cache_external_unpaired();

  ValueT  external_paired() const;
  ValueT& external_paired();
  void initialize_cache_external_paired();
#endif

private:
  static const std::string def_bases;
  static const size_t N;

  std::unordered_map<std::string,uint> hash_;
  std::vector<std::pair<std::string,ValueT>> keyval_;
  std::array<std::array<int, 256>, 256> is_complementary_;
  std::array<int, 256> is_base_;

  // cache
#ifdef USE_CACHE
#if PARAMS_BASE_PAIR
  VVI cache_base_pair_;
#endif
#if PARAMS_BASE_PAIR_DIST
  VI cache_base_pair_dist_at_least_;
#endif
#if PARAMS_TERMINAL_MISMATCH
  VVVVI cache_terminal_mismatch_;
#endif
#if PARAMS_HAIRPIN_LENGTH
  VI cache_hairpin_length_at_least_;
#endif
#if PARAMS_HELIX_LENGTH
  VI cache_helix_length_at_least_;
#endif
#if PARAMS_ISOLATED_BASE_PAIR
  int cache_isolated_base_pair_;
#endif
#if PARAMS_INTERNAL_EXPLICIT
  VVI cache_internal_explicit_;
#endif
#if PARAMS_BULGE_LENGTH
  VI cache_bulge_length_at_least_;
#endif
#if PARAMS_INTERNAL_LENGTH
  VI cache_internal_length_at_least_;
#endif
#if PARAMS_INTERNAL_SYMMETRY
  VI cache_internal_symmetric_length_at_least_;
#endif
#if PARAMS_INTERNAL_ASYMMETRY
  VI cache_internal_asymmetry_at_least_;
#endif
#if PARAMS_HELIX_STACKING
  VVVVI cache_helix_stacking_;
#endif
#if PARAMS_HELIX_CLOSING
  VVI cache_helix_closing_;
#endif
#if PARAMS_MULTI_LENGTH
  int cache_multi_base_;
  int cache_multi_unpaired_;
  int cache_multi_paired_;
#endif
#if PARAMS_DANGLE
  VVVI cache_dangle_left_;
  VVVI cache_dangle_right_;
#endif
#if PARAMS_EXTERNAL_LENGTH
  int cache_external_unpaired_;
  int cache_external_paired_;
#endif
#endif
};

#endif  // PARAMETERHASH_HPP

// Local Variables:
// mode: C++
// End:

