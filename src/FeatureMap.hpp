#ifndef FEATUREMAP_HPP
#define FEATUREMAP_HPP

#include "Config.hpp"
#include "Utilities.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <array>

#define USE_CACHE

class FeatureMap
{
public:
  FeatureMap(const char* def_bases = "ACGUP");
  ~FeatureMap() { }

  void initialize_cache();

  size_t set_key(const std::string& key);
  size_t find_key(const std::string& key) const
  {
    auto itr = hash_.find(key);
    return itr != hash_.end() ? itr->second : -1u;
  }

  // access to parameters
#if PARAMS_BASE_PAIR
  size_t base_pair(NUCL i, NUCL j) const;
  void initialize_cache_base_pair();
#endif
#if PARAMS_BASE_PAIR_DIST
  size_t  base_pair_dist_at_least(uint l) const;
  void initialize_cache_base_pair_dist_at_least();
#endif
#if PARAMS_TERMINAL_MISMATCH
  size_t terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const;
  void initialize_cache_terminal_mismatch();
#endif
#if PARAMS_HAIRPIN_LENGTH
  size_t hairpin_length_at_least(uint l) const;
  void initialize_cache_hairpin_length_at_least();
#endif

  size_t hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l) const;

#if PARAMS_HELIX_LENGTH
  size_t helix_length_at_least(uint l) const;
  void initialize_cache_helix_length_at_least();
#endif
#if PARAMS_ISOLATED_BASE_PAIR
  size_t isolated_base_pair() const;
  void initialize_cache_isolated_base_pair();
#endif
#if PARAMS_INTERNAL_EXPLICIT
  size_t internal_explicit(uint i, uint j) const;
  void initialize_cache_internal_explicit();
#endif
#if PARAMS_BULGE_LENGTH
  size_t bulge_length_at_least(uint l) const;
  void initialize_cache_bulge_length_at_least();
#endif
#if PARAMS_INTERNAL_LENGTH
  size_t internal_length_at_least(uint l) const;
  void initialize_cache_internal_length_at_least();
#endif
#if PARAMS_INTERNAL_SYMMETRY
  size_t internal_symmetric_length_at_least(uint l) const;
  void initialize_cache_internal_symmetric_length_at_least();
#endif
#if PARAMS_INTERNAL_ASYMMETRY
  size_t internal_asymmetry_at_least(uint l) const;
  void initialize_cache_internal_asymmetry_at_least();
#endif

  size_t internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m) const;

#if PARAMS_HELIX_STACKING
  size_t helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const;
  void initialize_cache_helix_stacking();
#endif
#if PARAMS_HELIX_CLOSING
  size_t helix_closing(NUCL i, NUCL j) const;
  void initialize_cache_helix_closing();
#endif
#if PARAMS_MULTI_LENGTH
  size_t multi_base() const;
  void initialize_cache_multi_base();

  size_t multi_unpaired() const;
  void initialize_cache_multi_unpaired();

  size_t multi_paired() const;
  void initialize_cache_multi_paired();
#endif
#if PARAMS_DANGLE
  size_t dangle_left(NUCL i1, NUCL j1, NUCL i2) const;
  void initialize_cache_dangle_left();

  size_t dangle_right(NUCL i1, NUCL j1, NUCL j2) const;
  void initialize_cache_dangle_right();
#endif
#if PARAMS_EXTERNAL_LENGTH
  size_t external_unpaired() const;
  void initialize_cache_external_unpaired();

  size_t external_paired() const;
  void initialize_cache_external_paired();
#endif

private:
  const std::string def_bases_;
  size_t NBASES;
  std::unordered_map<std::string, size_t> hash_;
  std::vector<const char*> keys_;
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

#endif  // FEATUREMAP_HPP
