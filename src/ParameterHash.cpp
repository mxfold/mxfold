#include "../config.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include <iterator>
#include <stdexcept>
#include <cerrno>
#include <cstdio>
#include <cctype>
#include <cassert>
#include "ParameterHash.hpp"
#include "Utilities.hpp"
#include "Config.hpp"

// static
template < class ValueT >
const std::string ParameterHash<ValueT>::def_bases("ACGUP");
template < class ValueT >
const size_t ParameterHash<ValueT>::N = 5; //ParameterHash<ValueT>::def_bases.size();

template < class ValueT >
ParameterHash<ValueT>::
ParameterHash()
  : trie_(), keys_(), values_()
#if PARAMS_BASE_PAIR
  , cache_base_pair_(N, VI(N, -1))
#endif
#if PARAMS_BASE_PAIR_DIST
  , cache_base_pair_dist_at_least_(D_MAX_BP_DIST_THRESHOLDS, -1)
#endif
#if PARAMS_TERMINAL_MISMATCH
  , cache_terminal_mismatch_(N, VVVI(N, VVI(N, VI(N, -1))))
#endif
#if PARAMS_HAIRPIN_LENGTH
  , cache_hairpin_length_at_least_(D_MAX_HAIRPIN_LENGTH, -1)
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
  , cache_hairpin_3_nucleotides_(N, VVI(N, VI(N, -1)))
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
  , cache_hairpin_4_nucleotides_(N, VVVI(N, VVI(N, VI(N, -1))))
#endif
#if PARAMS_HELIX_LENGTH
  , cache_helix_length_at_least_(D_MAX_HELIX_LENGTH, -1)
#endif
#if PARAMS_ISOLATED_BASE_PAIR
  , cache_isolated_base_pair_(-1)
#endif
#if PARAMS_INTERNAL_EXPLICIT
  , cache_internal_explicit_(D_MAX_INTERNAL_EXPLICIT_LENGTH+1, VI(D_MAX_INTERNAL_EXPLICIT_LENGTH+1, -1))
#endif
#if PARAMS_BULGE_LENGTH
  , cache_bulge_length_at_least_(D_MAX_BULGE_LENGTH, -1)
#endif
#if PARAMS_INTERNAL_LENGTH
  , cache_internal_length_at_least_(D_MAX_INTERNAL_LENGTH, -1)
#endif
#if PARAMS_INTERNAL_SYMMETRY
  , cache_internal_symmetric_length_at_least_(D_MAX_INTERNAL_SYMMETRIC_LENGTH, -1)
#endif
#if PARAMS_INTERNAL_ASYMMETRY
  , cache_internal_asymmetry_at_least_(D_MAX_INTERNAL_ASYMMETRY, -1)
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
  , cache_internal_0x1_nucleotides_(N, -1)
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
  , cache_internal_1x1_nucleotides_(N, VI(N, -1))
#endif
#if PARAMS_HELIX_STACKING
  , cache_helix_stacking_(N, VVVI(N, VVI(N, VI(N, -1))))
#endif
#if PARAMS_HELIX_CLOSING
  , cache_helix_closing_(N, VI(N, -1))
#endif
#if PARAMS_MULTI_LENGTH
  , cache_multi_base_(-1)
  , cache_multi_unpaired_(-1)
  , cache_multi_paired_(-1)
#endif
#if PARAMS_DANGLE
  , cache_dangle_left_(N, VVI(N, VI(N, -1)))
  , cache_dangle_right_(N, VVI(N, VI(N, -1)))
#endif
#if PARAMS_EXTERNAL_LENGTH
  , cache_external_unpaired_(-1)
  , cache_external_paired_(-1)
#endif
{
  initialize();
  initialize_cache();
}

template < class ValueT >
ParameterHash<ValueT>::
ParameterHash(ParameterHash&& other)
  : trie_(std::move(other.trie_)),
    keys_(std::move(other.keys_)),
    values_(std::move(other.values_))
#if PARAMS_BASE_PAIR
  , cache_base_pair_(std::move(other.cache_base_pair_))
#endif
#if PARAMS_BASE_PAIR_DIST
  , cache_base_pair_dist_at_least_(std::move(other.cache_base_pair_dist_at_least_))
#endif
#if PARAMS_TERMINAL_MISMATCH
  , cache_terminal_mismatch_(std::move(other.cache_terminal_mismatch_))
#endif
#if PARAMS_HAIRPIN_LENGTH
  , cache_hairpin_length_at_least_(std::move(other.cache_hairpin_length_at_least_))
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
  , cache_hairpin_3_nucleotides_(std::move(other.cache_hairpin_3_nucleotides_))
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
  , cache_hairpin_4_nucleotides_(std::move(other.cache_hairpin_4_nucleotides_))
#endif
#if PARAMS_HELIX_LENGTH
  , cache_helix_length_at_least_(std::move(other.cache_helix_length_at_least_))
#endif
#if PARAMS_ISOLATED_BASE_PAIR
  , cache_isolated_base_pair_(std::move(other.cache_isolated_base_pair_))
#endif
#if PARAMS_INTERNAL_EXPLICIT
  , cache_internal_explicit_(std::move(other.cache_internal_explicit_))
#endif
#if PARAMS_BULGE_LENGTH
  , cache_bulge_length_at_least_(std::move(other.cache_bulge_length_at_least_))
#endif
#if PARAMS_INTERNAL_LENGTH
  , cache_internal_length_at_least_(std::move(other.cache_internal_length_at_least_))
#endif
#if PARAMS_INTERNAL_SYMMETRY
  , cache_internal_symmetric_length_at_least_(std::move(other.cache_internal_symmetric_length_at_least_))
#endif
#if PARAMS_INTERNAL_ASYMMETRY
  , cache_internal_asymmetry_at_least_(std::move(other.cache_internal_asymmetry_at_least_))
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
  , cache_internal_0x1_nucleotides_(std::move(other.cache_internal_0x1_nucleotides_))
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
  , cache_internal_1x1_nucleotides_(std::move(other.cache_internal_1x1_nucleotides_))
#endif
#if PARAMS_HELIX_STACKING
  , cache_helix_stacking_(std::move(other.cache_helix_stacking_))
#endif
#if PARAMS_HELIX_CLOSING
  , cache_helix_closing_(std::move(other.cache_helix_closing_))
#endif
#if PARAMS_MULTI_LENGTH
  , cache_multi_base_(std::move(other.cache_multi_base_))
  , cache_multi_unpaired_(std::move(other.cache_multi_unpaired_))
  , cache_multi_paired_(std::move(other.cache_multi_paired_))
#endif
#if PARAMS_DANGLE
  , cache_dangle_left_(std::move(other.cache_dangle_left_))
  , cache_dangle_right_(std::move(other.cache_dangle_right_))
#endif
#if PARAMS_EXTERNAL_LENGTH
  , cache_external_unpaired_(std::move(other.cache_external_unpaired_))
  , cache_external_paired_(std::move(other.cache_external_paired_))
#endif
{
  initialize();
}


template < class ValueT >
ParameterHash<ValueT>&
ParameterHash<ValueT>::
operator=(ParameterHash&& other)
{
  trie_ = std::move(other.trie_);
  keys_ = std::move(other.keys_);
  values_ = std::move(other.values_);
#if PARAMS_BASE_PAIR
  cache_base_pair_ = std::move(other.cache_base_pair_);
#endif
#if PARAMS_BASE_PAIR_DIST
  cache_base_pair_dist_at_least_ = std::move(other.cache_base_pair_dist_at_least_);
#endif
#if PARAMS_TERMINAL_MISMATCH
  cache_terminal_mismatch_ = std::move(other.cache_terminal_mismatch_);
#endif
#if PARAMS_HAIRPIN_LENGTH
  cache_hairpin_length_at_least_ = std::move(other.cache_hairpin_length_at_least_);
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
  cache_hairpin_3_nucleotides_ = std::move(other.cache_hairpin_3_nucleotides_);
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
  cache_hairpin_4_nucleotides_ = std::move(other.cache_hairpin_4_nucleotides_);
#endif
#if PARAMS_HELIX_LENGTH
  cache_helix_length_at_least_ = std::move(other.cache_helix_length_at_least_);
#endif
#if PARAMS_ISOLATED_BASE_PAIR
  cache_isolated_base_pair_ = std::move(other.cache_isolated_base_pair_);
#endif
#if PARAMS_INTERNAL_EXPLICIT 
  cache_internal_explicit_ = std::move(other.cache_internal_explicit_);
#endif
#if PARAMS_BULGE_LENGTH
  cache_bulge_length_at_least_ = std::move(other.cache_bulge_length_at_least_);
#endif
#if PARAMS_INTERNAL_LENGTH
  cache_internal_length_at_least_ = std::move(other.cache_internal_length_at_least_);
#endif
#if PARAMS_INTERNAL_SYMMETRY
  cache_internal_symmetric_length_at_least_ = std::move(other.cache_internal_symmetric_length_at_least_);
#endif
#if PARAMS_INTERNAL_ASYMMETRY
  cache_internal_asymmetry_at_least_ = std::move(other.cache_internal_asymmetry_at_least_);
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
  cache_internal_0x1_nucleotides_ = std::move(other.cache_internal_0x1_nucleotides_);
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
  cache_internal_1x1_nucleotides_ = std::move(other.cache_internal_1x1_nucleotides_);
#endif
#if PARAMS_HELIX_STACKING
  cache_helix_stacking_ = std::move(other.cache_helix_stacking_);
#endif
#if PARAMS_HELIX_CLOSING
  cache_helix_closing_ = std::move(other.cache_helix_closing_);
#endif
#if PARAMS_MULTI_LENGTH
  cache_multi_base_ = std::move(other.cache_multi_base_);
  cache_multi_unpaired_ = std::move(other.cache_multi_unpaired_);
  cache_multi_paired_ = std::move(other.cache_multi_paired_);
#endif
#if PARAMS_DANGLE
  cache_dangle_left_ = std::move(other.cache_dangle_left_);
  cache_dangle_right_ = std::move(other.cache_dangle_right_);
#endif
#if PARAMS_EXTERNAL_LENGTH
  cache_external_unpaired_ = std::move(other.cache_external_unpaired_);
  cache_external_paired_ = std::move(other.cache_external_paired_);
#endif
  return *this;
}

template < class ValueT >
void
ParameterHash<ValueT>::
initialize()
{
  std::fill(std::begin(is_base_), std::end(is_base_), -1);
  for (size_t i=0; i!=def_bases.size(); ++i)
    is_base_[def_bases[i]] = i;

  // precompute complementary pairings
  for (auto e : is_complementary_)
    std::fill(std::begin(e), std::end(e), false);

  is_complementary_['A']['U'] = is_complementary_['U']['A'] = true;
  is_complementary_['G']['U'] = is_complementary_['U']['G'] = true;
  is_complementary_['C']['G'] = is_complementary_['G']['C'] = true;
  // Pseudouridine
  is_complementary_['A']['P'] = is_complementary_['P']['A'] = true;
  is_complementary_['C']['P'] = is_complementary_['P']['C'] = true;
  is_complementary_['G']['P'] = is_complementary_['P']['G'] = true;
  is_complementary_['U']['P'] = is_complementary_['P']['U'] = true;
}

template < class ValueT >
void
ParameterHash<ValueT>::
initialize_cache()
{
#if PARAMS_BASE_PAIR
  initialize_cache_base_pair();
#endif
#if PARAMS_BASE_PAIR_DIST
  initialize_cache_base_pair_dist_at_least();
#endif
#if PARAMS_TERMINAL_MISMATCH
  initialize_cache_terminal_mismatch();
#endif
#if PARAMS_HAIRPIN_LENGTH
  initialize_cache_hairpin_length_at_least();
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
  initialize_cache_hairpin_3_nucleotides();
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
  initialize_cache_hairpin_4_nucleotides();
#endif
#if PARAMS_HELIX_LENGTH
  initialize_cache_helix_length_at_least();
#endif
#if PARAMS_ISOLATED_BASE_PAIR
  initialize_cache_isolated_base_pair();
#endif
#if PARAMS_INTERNAL_EXPLICIT
  initialize_cache_internal_explicit();
#endif
#if PARAMS_BULGE_LENGTH
  initialize_cache_bulge_length_at_least();
#endif
#if PARAMS_INTERNAL_LENGTH
  initialize_cache_internal_length_at_least();
#endif
#if PARAMS_INTERNAL_SYMMETRY
  initialize_cache_internal_symmetric_length_at_least();
#endif
#if PARAMS_INTERNAL_ASYMMETRY
  initialize_cache_internal_asymmetry_at_least();
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
  initialize_cache_internal_0x1_nucleotides();
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
  initialize_cache_internal_1x1_nucleotides();
#endif
#if PARAMS_HELIX_STACKING
  initialize_cache_helix_stacking();
#endif
#if PARAMS_HELIX_CLOSING
  initialize_cache_helix_closing();
#endif
#if PARAMS_MULTI_LENGTH
  initialize_cache_multi_base();
  initialize_cache_multi_unpaired();
  initialize_cache_multi_paired();
#endif
#if PARAMS_DANGLE
  initialize_cache_dangle_left();
  initialize_cache_dangle_right();
#endif
#if PARAMS_EXTERNAL_LENGTH
  initialize_cache_external_unpaired();
  initialize_cache_external_paired();
#endif
}

template < class ValueT >
bool
ParameterHash<ValueT>::
is_complementary(NUCL x, NUCL y) const
{
  return is_complementary_[x][y];
}

template < class ValueT >
bool
ParameterHash<ValueT>::
is_base(NUCL x) const
{
  return is_base_[x]>=0;
}


template < class ValueT >
void
ParameterHash<ValueT>::
LoadFromHash(const std::unordered_map<std::string, ValueT>& hash)
{
  trie_.clear();
  keys_.clear();
  values_.clear();
  for (auto e: hash)
    set_key_value(e.first, e.second);
}

template < class ValueT >
void
ParameterHash<ValueT>::
ReadFromFile(const std::string& filename)
{
  trie_.clear();
  keys_.clear();
  values_.clear();
  std::ifstream is(filename.c_str());
  if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  std::string k;
  ValueT v;
  if (!std::isdigit(is.peek()))
  {
    while (is >> k >> v)
      if (v!=0.0)
        set_key_value(k, v);
  }
  else
  {
    // for reading AdaGradRDA outputs
    ValueT eta, lambda, eps, t, s1, s2;
    is >> eta >> lambda >> eps >> t;
    while (is >> k >> v >> s1 >> s2)
      if (v!=0.0)
        set_key_value(k, v);
  }
}

template < class ValueT >
void
ParameterHash<ValueT>::
WriteToFile(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
  if (!os) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  std::vector<size_t> idx(keys_.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&](size_t i, size_t j) { return keys_[i] < keys_[j]; });
  for (auto i: idx)
    if (values_[i]!=0.0)
      os << keys_[i] << " " << values_[i] << std::endl;
}

template < class ValueT >
inline
ValueT
ParameterHash<ValueT>::
get_by_key(const std::string& key) const
{
  auto i = trie_.exactMatchSearch<int>(key.c_str());
  return i==trie_t::CEDAR_NO_VALUE ? static_cast<ValueT>(0) : values_[i];
}

template < class ValueT >
inline
ValueT&
ParameterHash<ValueT>::
get_by_key(const std::string& key)
{
  auto i = trie_.exactMatchSearch<int>(key.c_str());
  if (i!=trie_t::CEDAR_NO_VALUE)
    return values_[i];

  trie_.update(key.c_str()) = values_.size();
  keys_.push_back(key);
  values_.push_back(static_cast<ValueT>(0.0));
  return values_.back();
}

template < class ValueT >
inline
uint
ParameterHash<ValueT>::
set_key_value(const std::string& key, ValueT val)
{
  auto i = set_key_value(key);
  values_[i] = val;
  return i;
}

template < class ValueT >
inline
uint
ParameterHash<ValueT>::
set_key_value(const std::string& key)
{
  auto i = trie_.exactMatchSearch<int>(key.c_str());
  if (i!=trie_t::CEDAR_NO_VALUE)
    return i;

  keys_.push_back(key);
  trie_.update(key.c_str()) = values_.size();

  // symmetric features
  const std::string s_internal_nucleotides("internal_nucleotides_");
  const std::string s_base_pair("base_pair_");
  const std::string s_helix_stacking("helix_stacking_");
  const std::string s_internal_explicit("internal_explicit_");
  if (key.find(s_internal_nucleotides) == 0)
  {
    size_t pos=key.find("_", s_internal_nucleotides.size());
    std::string nuc1 = key.substr(s_internal_nucleotides.size(), pos-s_internal_nucleotides.size());
    std::string nuc2 = key.substr(pos+1);
    std::reverse(nuc1.begin(), nuc1.end());
    std::reverse(nuc2.begin(), nuc2.end());
    std::string key2 = s_internal_nucleotides + nuc2 + "_" + nuc1;
    trie_.update(key2.c_str()) = values_.size();
    if (keys_.back()>key2) keys_.back() = key2;
  }
  else if (key.find(s_base_pair) == 0 && key.size() == s_base_pair.size()+2)
  {
    std::string nuc1 = key.substr(s_base_pair.size()+0, 1);
    std::string nuc2 = key.substr(s_base_pair.size()+1, 1);
    std::string key2 = s_base_pair + nuc2 + nuc1;
    trie_.update(key2.c_str()) = values_.size();
    if (keys_.back()>key2) keys_.back() = key2;
  }
  else if (key.find(s_helix_stacking) == 0)
  {
    std::string nuc1 = key.substr(s_helix_stacking.size()+0, 1);
    std::string nuc2 = key.substr(s_helix_stacking.size()+1, 1);
    std::string nuc3 = key.substr(s_helix_stacking.size()+2, 1);
    std::string nuc4 = key.substr(s_helix_stacking.size()+3, 1);
    std::string key2 = s_helix_stacking + nuc4 + nuc3 + nuc2 + nuc1;
    trie_.update(key2.c_str()) = values_.size();
    if (keys_.back()>key2) keys_.back() = key2;
  }
  else if (key.find(s_internal_explicit) == 0)
  {
    size_t pos=key.find("_", s_internal_explicit.size());
    std::string l1 = key.substr(s_internal_explicit.size(), pos-s_internal_explicit.size());
    std::string l2 = key.substr(pos+1);
    std::string key2 = s_internal_explicit + l2 + "_" + l1;
    trie_.update(key2.c_str()) = values_.size();
    if (keys_.back()>key2) keys_.back() = key2;
  }

  values_.push_back(static_cast<ValueT>(0));
  return values_.size()-1;
}

#if PARAMS_BASE_PAIR
static const char* format_base_pair = "base_pair_%c%c";

template < class ValueT >
inline
void
ParameterHash<ValueT>::
initialize_cache_base_pair()
{
  for (auto i : def_bases)
    for (auto j : def_bases)
      cache_base_pair_[is_base_[i]][is_base_[j]] = set_key_value(SPrintF(format_base_pair, i, j));
}

template < class ValueT >
inline
ValueT
ParameterHash<ValueT>::
base_pair(NUCL i, NUCL j) const
{
  auto ii=is_base_[i], jj=is_base_[j];
  if (ii>=0 && jj>=0)
    return values_[cache_base_pair_[ii][jj]];
  else
    return get_by_key(SPrintF(format_base_pair, i, j));
}

template < class ValueT >
inline
ValueT&
ParameterHash<ValueT>::
base_pair(NUCL i, NUCL j)
{
  auto ii=is_base_[i], jj=is_base_[j];
  if (ii>=0 && jj>=0)
    return values_[cache_base_pair_[ii][jj]];
  else
    return get_by_key(SPrintF(format_base_pair, i, j));
}
#endif

#if PARAMS_BASE_PAIR_DIST
static const char* format_base_pair_dist_at_least = "base_pair_dist_at_least_%d";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_base_pair_dist_at_least()
{
  for (size_t i=0; i!=cache_base_pair_dist_at_least_.size(); ++i)
    cache_base_pair_dist_at_least_[i] = set_key_value(SPrintF(format_base_pair_dist_at_least, i));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
base_pair_dist_at_least(uint l) const
{
  if (l<cache_base_pair_dist_at_least_.size())
    return values_[cache_base_pair_dist_at_least_[l]];
  else
    return get_by_key(SPrintF(format_base_pair_dist_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
base_pair_dist_at_least(uint l)
{
  if (l<cache_base_pair_dist_at_least_.size())
    return values_[cache_base_pair_dist_at_least_[l]];
  else
    return get_by_key(SPrintF(format_base_pair_dist_at_least, l));
}
#endif

#if PARAMS_TERMINAL_MISMATCH
static const char* format_terminal_mismatch = "terminal_mismatch_%c%c%c%c";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_terminal_mismatch()
{
  for (auto i1: def_bases)
  {
    auto ii1 = is_base_[i1];
    for (auto j1: def_bases)
    {
      auto jj1 = is_base_[j1];
      for (auto i2: def_bases)
      {
        auto ii2 = is_base_[i2];
        for (auto j2: def_bases)
        {
          auto jj2 = is_base_[j2];
          cache_terminal_mismatch_[ii1][jj1][ii2][jj2]
            = set_key_value(SPrintF(format_terminal_mismatch, i1, j1, i2, j2));
        }
      }
    }
  }
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
  auto ii1 = is_base_[i1], jj1 = is_base_[j1];
  auto ii2 = is_base_[i2], jj2 = is_base_[j2];
  if (ii1>=0 && jj1>=0 && ii2>=0 && jj2>=0)
    return values_[cache_terminal_mismatch_[ii1][jj1][ii2][jj2]];
  else
    return get_by_key(SPrintF(format_terminal_mismatch, i1, j1, i2, j2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2)
{
  auto ii1 = is_base_[i1], jj1 = is_base_[j1];
  auto ii2 = is_base_[i2], jj2 = is_base_[j2];
  if (ii1>=0 && jj1>=0 && ii2>=0 && jj2>=0)
    return values_[cache_terminal_mismatch_[ii1][jj1][ii2][jj2]];
  else
    return get_by_key(SPrintF(format_terminal_mismatch, i1, j1, i2, j2));
}
#endif

#if PARAMS_HAIRPIN_LENGTH
static const char* format_hairpin_length_at_least = "hairpin_length_at_least_%d";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_hairpin_length_at_least()
{
  for (size_t i=0; i!=cache_hairpin_length_at_least_.size(); ++i)
    cache_hairpin_length_at_least_[i] = set_key_value(SPrintF(format_hairpin_length_at_least, i));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_length_at_least(uint l) const
{
  if (l<cache_hairpin_length_at_least_.size())
    return values_[cache_hairpin_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_hairpin_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_length_at_least(uint l)
{
  if (l<cache_hairpin_length_at_least_.size())
    return values_[cache_hairpin_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_hairpin_length_at_least, l));
}
#endif

static const std::string format_hairpin_nucleotides("hairpin_nucleotides_");

#if PARAMS_HAIRPIN_3_NUCLEOTIDES
template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_hairpin_3_nucleotides()
{
  for (auto i: def_bases)
  {
    auto ii = is_base_[i];
    for (auto j: def_bases)
    {
      auto jj = is_base_[j];
      for (auto k: def_bases)
      {
        auto kk = is_base_[k];
        std::string key3 = format_hairpin_nucleotides+i+j+k;
        cache_hairpin_3_nucleotides_[ii][jj][kk] = set_key_value(key3);
      }
    }
  }
}
#endif

#if PARAMS_HAIRPIN_4_NUCLEOTIDES
template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_hairpin_4_nucleotides()
{
  for (auto i: def_bases)
  {
    auto ii = is_base_[i];
    for (auto j: def_bases)
    {
      auto jj = is_base_[j];
      for (auto k: def_bases)
      {
        auto kk = is_base_[k];
        std::string key3 = format_hairpin_nucleotides+i+j+k;
        for (auto l: def_bases)
        {
          auto ll = is_base_[l];
          cache_hairpin_4_nucleotides_[ii][jj][kk][ll] = set_key_value(key3+l);
        }
      }
    }
  }
}
#endif

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l) const
{
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
  if (l==3)
  {
    auto ii = is_base_[s[i+0]];
    auto jj = is_base_[s[i+1]];
    auto kk = is_base_[s[i+2]];
    return values_[cache_hairpin_3_nucleotides_[ii][jj][kk]];
  }
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
  if (l==4)
  {
    auto ii = is_base_[s[i+0]];
    auto jj = is_base_[s[i+1]];
    auto kk = is_base_[s[i+2]];
    auto ll = is_base_[s[i+3]];
    return values_[cache_hairpin_4_nucleotides_[ii][jj][kk][ll]];
  }
#endif

  std::string h(l, ' ');
  std::copy(&s[i], &s[i]+l, h.begin());
  return get_by_key(format_hairpin_nucleotides + h);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l)
{
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
  if (l==3)
  {
    auto ii = is_base_[s[i+0]];
    auto jj = is_base_[s[i+1]];
    auto kk = is_base_[s[i+2]];
    return values_[cache_hairpin_3_nucleotides_[ii][jj][kk]];
  }
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
  if (l==4)
  {
    auto ii = is_base_[s[i+0]];
    auto jj = is_base_[s[i+1]];
    auto kk = is_base_[s[i+2]];
    auto ll = is_base_[s[i+3]];
    return values_[cache_hairpin_4_nucleotides_[ii][jj][kk][ll]];
  }
#endif

  std::string h(l, ' ');
  auto x = std::copy(&s[i], &s[i+l], h.begin());
  return get_by_key(format_hairpin_nucleotides + h);
}

template <class ValueT>
inline
VI
ParameterHash<ValueT>::
hairpin_nucleotides_cache(const std::vector<NUCL>& s, uint i, uint max_l) const
{
  static size_t head = -1u;
  VI v(max_l, trie_t::CEDAR_NO_VALUE);
  size_t node_pos=0, key_pos=0;
  if (head == -1u)
  {
    auto r = trie_.traverse(format_hairpin_nucleotides.c_str(), node_pos, key_pos, 
                            format_hairpin_nucleotides.size());
    if (r == trie_t::CEDAR_NO_PATH) return v;
    head = node_pos;
  }

  node_pos = head;
  for (size_t j=0; j!=max_l; ++j)
  {
    v[j] = trie_.traverse(&s[i], node_pos, key_pos, j+1);
    if (v[j]==trie_t::CEDAR_NO_PATH) break;
  }
  return v;
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l, const VI& pos) const
{
  return pos[l]>=0 ? values_[pos[l]] : static_cast<ValueT>(0);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l, const VI& pos)
{
#ifdef TEST_CACHE
  auto v = pos[l]>=0 ? values_[pos[l]] : hairpin_nucleotides(s, i, l);
  assert(v == hairpin_nucleotides(s, i, l));
#endif
  return pos[l]>=0 ? values_[pos[l]] : hairpin_nucleotides(s, i, l);
}

#if PARAMS_HELIX_LENGTH
static const char* format_helix_length_at_least = "helix_length_at_least_%d";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_helix_length_at_least()
{
  for (size_t i=0; i!=cache_helix_length_at_least_.size(); ++i)
    cache_helix_length_at_least_[i] = set_key_value(SPrintF(format_helix_length_at_least, i));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
helix_length_at_least(uint l) const
{
  if (l<cache_helix_length_at_least_.size())
    return values_[cache_helix_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_helix_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
helix_length_at_least(uint l)
{
  if (l<cache_helix_length_at_least_.size())
    return values_[cache_helix_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_helix_length_at_least, l));
}
#endif

#if PARAMS_ISOLATED_BASE_PAIR
static const char* format_isolated_base_pair = "isolated_base_pair";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_isolated_base_pair()
{
  cache_isolated_base_pair_ = set_key_value(format_isolated_base_pair);
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
isolated_base_pair() const
{
  return values_[cache_isolated_base_pair_];
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
isolated_base_pair()
{
  return values_[cache_isolated_base_pair_];
}
#endif

#if PARAMS_INTERNAL_EXPLICIT
static const char* format_internal_explicit = "internal_explicit_%d_%d";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_internal_explicit()
{
  for (size_t i=0; i!=cache_internal_explicit_.size(); ++i)
    for (size_t j=0; j!=cache_internal_explicit_[i].size(); ++j)
      cache_internal_explicit_[i][j] = set_key_value(SPrintF(format_internal_explicit, i, j));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_explicit(uint i, uint j) const
{
  if (i<cache_internal_explicit_.size() && j<cache_internal_explicit_[i].size())
    return values_[cache_internal_explicit_[i][j]];
  else
    return get_by_key(SPrintF(format_internal_explicit, i, j));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_explicit(uint i, uint j)
{
  if (i<cache_internal_explicit_.size() && j<cache_internal_explicit_[i].size())
    return values_[cache_internal_explicit_[i][j]];
  else
    return get_by_key(SPrintF(format_internal_explicit, i, j));
}
#endif

#if PARAMS_BULGE_LENGTH
static const char* format_bulge_length_at_least = "bulge_length_at_least_%d";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_bulge_length_at_least()
{
  for (size_t i=0; i!=cache_bulge_length_at_least_.size(); ++i)
    cache_bulge_length_at_least_[i] = set_key_value(SPrintF(format_bulge_length_at_least, i));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_length_at_least(uint l) const
{
  if (l<cache_bulge_length_at_least_.size())
    return values_[cache_bulge_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_bulge_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_length_at_least(uint l)
{
  if (l<cache_bulge_length_at_least_.size())
    return values_[cache_bulge_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_bulge_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_LENGTH
static const char* format_internal_length_at_least = "internal_length_at_least_%d";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_internal_length_at_least()
{
  for (size_t i=0; i!=cache_internal_length_at_least_.size(); ++i)
    cache_internal_length_at_least_[i] = set_key_value(SPrintF(format_internal_length_at_least, i));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_length_at_least(uint l) const
{
  if (l<cache_internal_length_at_least_.size())
    return values_[cache_internal_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_internal_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_length_at_least(uint l)
{
  if (l<cache_internal_length_at_least_.size())
    return values_[cache_internal_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_internal_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_SYMMETRY
static const char* format_internal_symmetric_length_at_least = "internal_symmetric_length_at_least_%d";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_internal_symmetric_length_at_least()
{
  for (size_t i=0; i!=cache_internal_symmetric_length_at_least_.size(); ++i)
    cache_internal_symmetric_length_at_least_[i] = set_key_value(SPrintF(format_internal_symmetric_length_at_least, i));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_symmetric_length_at_least(uint l) const
{
  if (l<cache_internal_symmetric_length_at_least_.size())
    return values_[cache_internal_symmetric_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_internal_symmetric_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_symmetric_length_at_least(uint l)
{
  if (l<cache_internal_symmetric_length_at_least_.size())
    return values_[cache_internal_symmetric_length_at_least_[l]];
  else
    return get_by_key(SPrintF(format_internal_symmetric_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_ASYMMETRY
static const char* format_internal_asymmetry_at_least = "internal_asymmetry_at_least_%d";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_internal_asymmetry_at_least()
{
  for (size_t i=0; i!=cache_internal_asymmetry_at_least_.size(); ++i)
    cache_internal_asymmetry_at_least_[i] = set_key_value(SPrintF(format_internal_asymmetry_at_least, i));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_asymmetry_at_least(uint l) const
{
  if (l<cache_internal_asymmetry_at_least_.size())
    return values_[cache_internal_asymmetry_at_least_[l]];
  else
    return get_by_key(SPrintF(format_internal_asymmetry_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_asymmetry_at_least(uint l)
{
  if (l<cache_internal_asymmetry_at_least_.size())
    return values_[cache_internal_asymmetry_at_least_[l]];
  else
    return get_by_key(SPrintF(format_internal_asymmetry_at_least, l));
}
#endif

static const std::string format_internal_nucleotides("internal_nucleotides_");

#if PARAMS_BULGE_0x1_NUCLEOTIDES
template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_internal_0x1_nucleotides()
{
  for (auto i1: def_bases)
  {
    auto ii1 = is_base_[i1];
    cache_internal_0x1_nucleotides_[ii1] = set_key_value(format_internal_nucleotides+"_"+i1);
  }
}
#endif

#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_internal_1x1_nucleotides()
{
  for (auto i1: def_bases)
  {
    auto ii1 = is_base_[i1];
    std::string key = format_internal_nucleotides+i1+"_";
    for (auto j1: def_bases)
    {
      auto jj1 = is_base_[j1];
      cache_internal_1x1_nucleotides_[ii1][jj1] = set_key_value(key+j1);
    }
  }
}
#endif

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m) const
{
#if PARAMS_BULGE_0x1_NUCLEOTIDES
  if (l==1 && m==0)
  {
    auto ii1 = is_base_[s[i]];
    return values_[cache_internal_0x1_nucleotides_[ii1]];
  }
  else if (l==0 && m==1)
  {
    auto jj1 = is_base_[s[j]];
    return values_[cache_internal_0x1_nucleotides_[jj1]];
  }
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
  if (l==1 && m==1)
  {
    auto ii1 = is_base_[s[i]];
    auto jj1 = is_base_[s[j]];
    return values_[cache_internal_1x1_nucleotides_[ii1][jj1]];
  }
#endif

  std::string nuc(l+m+1, ' ');
  auto x = std::copy(&s[i], &s[i+l], nuc.begin());
  *(x++) = '_';
  std::reverse_copy(&s[j-m+1], &s[j+1], x);
  return get_by_key(format_internal_nucleotides + nuc);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m)
{
#if PARAMS_BULGE_0x1_NUCLEOTIDES
  if (l==1 && m==0)
  {
    auto ii1 = is_base_[s[i]];
    return values_[cache_internal_0x1_nucleotides_[ii1]];
  }
  else if (l==0 && m==1)
  {
    auto jj1 = is_base_[s[j]];
    return values_[cache_internal_0x1_nucleotides_[jj1]];
  }
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
  if (l==1 && m==1)
  {
    auto ii1 = is_base_[s[i]];
    auto jj1 = is_base_[s[j]];
    return values_[cache_internal_1x1_nucleotides_[ii1][jj1]];
  }
#endif

  std::string nuc1(l+m+1, ' ');
  auto x = std::copy(&s[i], &s[i+l], nuc1.begin());
  *(x++) = '_';
  std::reverse_copy(&s[j-m+1], &s[j+1], x);

  auto key1 = format_internal_nucleotides+nuc1;
  auto k = trie_.exactMatchSearch<int>(key1.c_str());
  if (k!=trie_t::CEDAR_NO_VALUE)
    return values_[k];

  keys_.push_back(key1);
  trie_.update(key1.c_str()) = values_.size();

  std::string nuc2(l+m+1, ' ');
  x = std::copy(&s[j-m+1], &s[j+1], nuc2.begin());
  *(x++) = '_';
  std::reverse_copy(&s[i], &s[i+l], x);
  auto key2 = format_internal_nucleotides+nuc2;
  trie_.update(key2.c_str()) = values_.size();
  if (keys_.back()>key2) keys_.back() = key2;

  values_.push_back(static_cast<ValueT>(0.0));
  return values_.back();
}

template <class ValueT>
VVI
ParameterHash<ValueT>::
internal_nucleotides_cache(const std::vector<NUCL>& s, uint i, uint j,
                           uint max_l, uint max_m) const
{
  static size_t head = -1u;

  VVI ret(max_l+1, VI(max_m+1, trie_t::CEDAR_NO_VALUE));

  size_t node_pos=0, key_pos=0;
  int k;
  if (head == -1u)
  {
    k = trie_.traverse(format_internal_nucleotides.c_str(), node_pos, key_pos, 
                       format_internal_nucleotides.size());
    head = node_pos;
  }
  node_pos = head;

  for (uint l=0; l<=max_l && i+l-1<s.size(); ++l)
  {
    if (l>0)
    {
      key_pos = 0;
      k = trie_.traverse(&s[i+l-1], node_pos, key_pos, 1);
#ifdef TEST_CACHE
      std::string s2 = format_internal_nucleotides;
      for (uint kk=i; kk<i+l; ++kk) s2.push_back(s[kk]);
      size_t np=0, kp=0;
      auto k2 = trie_.traverse(s2.c_str(), np, kp, s2.size());
      assert(k == k2);
#endif
      if (k==trie_t::CEDAR_NO_PATH) break;
    }

    auto node_pos2 = node_pos;
    key_pos = 0;
    k = trie_.traverse("_", node_pos2, key_pos, 1);

    for (uint m=0; m<=max_m && j+1>=m && l+m<=DEFAULT_C_MAX_SINGLE_LENGTH; ++m)
    {
      if (l+m<1) continue;
      if (m>0)
      {
        key_pos = 0;
        k = trie_.traverse(&s[j-m+1], node_pos2, key_pos, 1);
#ifdef TEST_CACHE
        std::string s2 = format_internal_nucleotides;
        for (uint kk=i; kk<i+l; ++kk) s2.push_back(s[kk]);
        s2.push_back('_');
        for (uint kk=j; j+1<=kk+m && kk!=-1u; --kk) s2.push_back(s[kk]);
        size_t np=0, kp=0;
        auto k2 = trie_.traverse(s2.c_str(), np, kp, s2.size());
        assert(k == k2);
#endif
      }
      if (k==trie_t::CEDAR_NO_PATH) break;
      ret[l][m] = k;
    }
  }
  return ret;
}


template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m, const VVI& pos) const
{
#ifdef TEST_CACHE
  ValueT v;
  if (l<pos.size() && m<pos[l].size())
    v = pos[l][m]>=0 ? values_[pos[l][m]] : static_cast<ValueT>(0);
  else
    v = internal_nucleotides(s, i, l, j, m);
  assert(v == internal_nucleotides(s, i, l, j, m));
#endif
  if (l<pos.size() && m<pos[l].size())
        return pos[l][m]>=0 ? values_[pos[l][m]] : static_cast<ValueT>(0);
  else
    return internal_nucleotides(s, i, l, j, m);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m, const VVI& pos)
{
  if (l<pos.size() && m<pos[l].size())
    return pos[l][m]>=0 ? values_[pos[l][m]] : internal_nucleotides(s, i, l, j, m);
  else
    return internal_nucleotides(s, i, l, j, m);
}


#if PARAMS_HELIX_STACKING
static const char* format_helix_stacking = "helix_stacking_%c%c%c%c";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_helix_stacking()
{
  for (auto i1: def_bases)
  {
    auto ii1 = is_base_[i1];
    for (auto j1: def_bases)
    {
      auto jj1 = is_base_[j1];
      for (auto i2: def_bases)
      {
        auto ii2 = is_base_[i2];
        for (auto j2: def_bases)
        {
          auto jj2 = is_base_[j2];
          cache_helix_stacking_[ii1][jj1][ii2][jj2] = set_key_value(SPrintF(format_helix_stacking, i1, j1, i2, j2));
        }
      }
    }
  }
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
  auto ii1=is_base_[i1], jj1=is_base_[j1];
  auto ii2=is_base_[i2], jj2=is_base_[j2];
  if (ii1>=0 && jj1>=0 && ii2>=0 && jj2>=0)
    return values_[cache_helix_stacking_[ii1][jj1][ii2][jj2]];
  else
    return get_by_key(SPrintF(format_helix_stacking, i1, j1, i2, j2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2)
{
  auto ii1=is_base_[i1], jj1=is_base_[j1];
  auto ii2=is_base_[i2], jj2=is_base_[j2];
  if (ii1>=0 && jj1>=0 && ii2>=0 && jj2>=0)
    return values_[cache_helix_stacking_[ii1][jj1][ii2][jj2]];
  else
    return get_by_key(SPrintF(format_helix_stacking, i1, j1, i2, j2));
}
#endif

#if PARAMS_HELIX_CLOSING
static const char* format_helix_closing = "helix_closing_%c%c";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_helix_closing()
{
  for (auto i: def_bases)
  {
    auto ii = is_base_[i];
    for (auto j: def_bases)
    {
      auto jj = is_base_[j];
      cache_helix_closing_[ii][jj] = set_key_value(SPrintF(format_helix_closing, i, j));
    }
  }
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
helix_closing(NUCL i, NUCL j) const
{
  auto ii=is_base_[i], jj=is_base_[j];
  if (ii>=0 && jj>=0)
    return values_[cache_helix_closing_[ii][jj]];
  else
    return get_by_key(SPrintF(format_helix_closing, i, j));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
helix_closing(NUCL i, NUCL j)
{
  auto ii=is_base_[i], jj=is_base_[j];
  if (ii>=0 && jj>=0)
    return values_[cache_helix_closing_[ii][jj]];
  else
    return get_by_key(SPrintF(format_helix_closing, i, j));
}
#endif

#if PARAMS_MULTI_LENGTH
static const char* format_multi_base = "multi_base";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_multi_base()
{
  cache_multi_base_ = set_key_value(format_multi_base);
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
multi_base() const
{
  return values_[cache_multi_base_];
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
multi_base()
{
  return values_[cache_multi_base_];
}

static const char* format_multi_unpaired = "multi_unpaired";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_multi_unpaired()
{
  cache_multi_unpaired_ = set_key_value(format_multi_unpaired);
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
multi_unpaired() const
{
  return values_[cache_multi_unpaired_];
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
multi_unpaired()
{
  return values_[cache_multi_unpaired_];
}

static const char* format_multi_paired = "multi_paired";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_multi_paired()
{
  cache_multi_paired_ = set_key_value(format_multi_paired);
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
multi_paired() const
{
  return values_[cache_multi_paired_];
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
multi_paired()
{
  return values_[cache_multi_paired_];
}
#endif

#if PARAMS_DANGLE
static const char* format_dangle_left = "dangle_left_%c%c%c";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_dangle_left()
{
  for (auto i1: def_bases)
  {
    auto ii1 = is_base_[i1];
    for (auto j1: def_bases)
    {
      auto jj1 = is_base_[j1];
      for (auto i2: def_bases)
      {
        auto ii2 = is_base_[i2];
        cache_dangle_left_[ii1][jj1][ii2] = set_key_value(SPrintF(format_dangle_left, i1, j1, i2));
      }
    }
  }
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
dangle_left(NUCL i1, NUCL j1, NUCL i2) const
{
  auto ii1=is_base_[i1], jj1=is_base_[j1], ii2=is_base_[i2];
  if (ii1>=0 && jj1>=0 && ii2>=0)
    return values_[cache_dangle_left_[ii1][jj1][ii2]];
  else
    return get_by_key(SPrintF(format_dangle_left, i1, j1, i2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
dangle_left(NUCL i1, NUCL j1, NUCL i2)
{
  auto ii1=is_base_[i1], jj1=is_base_[j1], ii2=is_base_[i2];
  if (ii1>=0 && jj1>=0 && ii2>=0)
    return values_[cache_dangle_left_[ii1][jj1][ii2]];
  else
    return get_by_key(SPrintF(format_dangle_left, i1, j1, i2));
}

static const char* format_dangle_right = "dangle_right_%c%c%c";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_dangle_right()
{
  for (auto i1: def_bases)
  {
    auto ii1 = is_base_[i1];
    for (auto j1: def_bases)
    {
      auto jj1 = is_base_[j1];
      for (auto j2: def_bases)
      {
        auto jj2 = is_base_[j2];
        cache_dangle_right_[ii1][jj1][jj2] = set_key_value(SPrintF(format_dangle_right, i1, j1, j2));
      }
    }
  }
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
dangle_right(NUCL i1, NUCL j1, NUCL j2) const
{
  auto ii1=is_base_[i1], jj1=is_base_[j1], jj2=is_base_[j2];
  if (ii1>=0 && jj1>=0 && jj2>=0)
    return values_[cache_dangle_right_[ii1][jj1][jj2]];
  else
    return get_by_key(SPrintF(format_dangle_right, i1, j1, j2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
dangle_right(NUCL i1, NUCL j1, NUCL j2)
{
  auto ii1=is_base_[i1], jj1=is_base_[j1], jj2=is_base_[j2];
  if (ii1>=0 && jj1>=0 && jj2>=0)
    return values_[cache_dangle_right_[ii1][jj1][jj2]];
  else
    return get_by_key(SPrintF(format_dangle_right, i1, j1, j2));
}
#endif

#if PARAMS_EXTERNAL_LENGTH
static const char* format_external_unpaired = "external_unpaired";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_external_unpaired()
{
  cache_external_unpaired_ = set_key_value(format_external_unpaired);
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
external_unpaired() const
{
  return values_[cache_external_unpaired_];
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
external_unpaired()
{
  return values_[cache_external_unpaired_];
}

static const char* format_external_paired = "external_paired";

template <class ValueT>
inline
void
ParameterHash<ValueT>::
initialize_cache_external_paired()
{
  cache_external_paired_ = set_key_value(format_external_paired);
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
external_paired() const
{
  return values_[cache_external_paired_];
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
external_paired()
{
  return values_[cache_external_paired_];
}
#endif

template <class ValueT>
inline
bool
ParameterHash<ValueT>::
is_basepair_feature(const std::string& f) const
{
  static const char* bp_features[] = {
#if PARAMS_BASE_PAIR
    format_base_pair,
#endif
#if PARAMS_BASE_PAIR_DIST
    format_base_pair_dist_at_least,
#endif
#if PARAMS_TERMINAL_MISMATCH
    format_terminal_mismatch,
#endif
#if PARAMS_HELIX_LENGTH
    format_helix_length_at_least,
#endif
#if PARAMS_ISOLATED_BASE_PAIR
    format_isolated_base_pair,
#endif
#if PARAMS_HELIX_STACKING
    format_helix_stacking,
#endif
#if PARAMS_HELIX_CLOSING
    format_helix_closing,
#endif
#if PARAMS_DANGLE
    format_dangle_left, 
    format_dangle_right,
#endif
#if PARAMS_MULTI_LENGTH
    format_multi_paired,
#endif
#if PARAMS_EXTERNAL_LENGTH
    format_external_paired,
#endif
    NULL
  };

  for (uint i=0; bp_features[i]!=NULL; ++i)
  {
    auto m = std::mismatch(f.begin(), f.end(), bp_features[i]);
    if (*m.second == 0 || *m.second == '%')
      return true;
  }
  return false;
}

template <class ValueT>
inline
bool
ParameterHash<ValueT>::
is_basepair_context_feature(const std::string& f) const
{
  static const char* bp_features[] = {
#if PARAMS_BASE_PAIR
    format_base_pair,
#endif  
#if PARAMS_TERMINAL_MISMATCH
    format_terminal_mismatch,
#endif
#if PARAMS_HELIX_STACKING
    format_helix_stacking,
#endif
#if PARAMS_HELIX_CLOSING
    format_helix_closing,
#endif
#if PARAMS_DANGLE
    format_dangle_left, 
    format_dangle_right,
#endif
    NULL
  };

  for (uint i=0; bp_features[i]!=NULL; ++i)
  {
    auto m = std::mismatch(f.begin(), f.end(), bp_features[i]);
    if (is_base(*m.first) && *m.second == '%')
      return true;
  }
  return false;
}

template <class ValueT>
inline
bool
ParameterHash<ValueT>::
is_context_feature(const std::string& f) const
{
  static const char* context_features[] = {
    "hairpin_", "bulge_", "internal_", NULL
  };

  for (uint i=0; context_features[i]!=NULL; ++i)
  {
    auto m = std::mismatch(f.begin(), f.end(), context_features[i]);
    if (isdigit(*m.first))
      return true;
  }
  return false;
}

template
class ParameterHash<param_value_type>;

//template
//class ParameterHash<uint>;
