#include "../config.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include <regex>
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
std::string ParameterHash<ValueT>::def_bases = "ACGUP";


template<typename ... Args>
std::string
string_format( const std::string& format, Args ... args )
{
  size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
  std::unique_ptr<char[]> buf( new char[ size ] ); 
  snprintf( buf.get(), size, format.c_str(), args ... );
  return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
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
  for (auto e: hash)
  {
    set_key_value(e.first, e.second);
  }
}

template < class ValueT >
void
ParameterHash<ValueT>::
ReadFromFile(const std::string& filename)
{
  trie_.clear();
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

  std::vector<std::string> keys(values_.size());
  for (auto p=this->cbegin(); p!=this->cend(); ++p)
    if (keys[p.index()].empty() || p.key()<keys[p.index()])
        keys[p.index()] = p.key();
  std::vector<size_t> idx(keys.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), 
            [&keys](size_t i, size_t j) { return keys[i] < keys[j]; });
  for (auto i: idx)
    if (values_[i]!=0.0)
      os << keys[i] << " " << values_[i] << std::endl;
}

template < class ValueT >
inline
ValueT
ParameterHash<ValueT>::
get_by_key(const std::string& key) const
{
  auto i = trie_.template exactMatchSearch<int>(key.c_str());
  return i==trie_t::CEDAR_NO_VALUE ? static_cast<ValueT>(0) : values_[i];
}

template < class ValueT >
inline
ValueT&
ParameterHash<ValueT>::
get_by_key(const std::string& key)
{
  auto i = trie_.template exactMatchSearch<int>(key.c_str());
  if (i!=trie_t::CEDAR_NO_VALUE)
    return values_[i];

  trie_.update(key.c_str()) = values_.size();
  values_.push_back(static_cast<ValueT>(0.0));
  return values_.back();
}

template < class ValueT >
inline
void
ParameterHash<ValueT>::
set_key_value(const std::string& key, ValueT val)
{
  auto i = trie_.template exactMatchSearch<int>(key.c_str());
  if (i!=trie_t::CEDAR_NO_VALUE)
  {
    values_[i] = val;
    return;
  }

  trie_.update(key.c_str()) = values_.size();

  // symmetric features
  std::regex re0("internal_nucleotides_([A-Za-z]*)_([A-Za-z]*)");
  std::regex re1("base_pair_([A-Za-z])([A-Za-z])");
  std::regex re2("helix_stacking_([A-Za-z])([A-Za-z])([A-Za-z])([A-Za-z])");
  std::regex re3("internal_explicit_([0-9]+)_([0-9]+)");
  std::smatch match;
  if (std::regex_match(key, match, re0)) // internal_nucleotides
  {
    std::string nuc1(match[1]), nuc2(match[2]);
    std::reverse(nuc1.begin(), nuc1.end());
    std::reverse(nuc2.begin(), nuc2.end());
    std::string key2("internal_nucleotides_");
    key2 += nuc2 + "_" + nuc1;
    trie_.update(key2.c_str()) = values_.size();
  }
  else if (std::regex_match(key, match, re1)) // base_pair
  {
    std::string nuc1(match[1]), nuc2(match[2]);
    std::string key2("base_pair_");
    key2 += nuc2 + nuc1;
    trie_.update(key2.c_str()) = values_.size();
  }
  else if (std::regex_match(key, match, re2)) // helix_stacking
  {
    std::string nuc1(match[1]), nuc2(match[2]), nuc3(match[3]), nuc4(match[4]);
    std::string key2("helix_stacking_");
    key2 += nuc4 + nuc3 + nuc2 + nuc1;
    trie_.update(key2.c_str()) = values_.size();
  }
  else if (std::regex_match(key, match, re3)) // internal_explicit
  {
    std::string l1(match[1]), l2(match[2]);
    std::string key2("internal_explicit_");
    key2 += l2 + "_" + l1;
    trie_.update(key2.c_str()) = values_.size();
  }

  values_.push_back(val);
}

#if PARAMS_BASE_PAIR
static const char* format_base_pair = "base_pair_%c%c";

template < class ValueT >
inline
ValueT
ParameterHash<ValueT>::
base_pair(NUCL i, NUCL j) const
{
  return get_by_key(string_format(format_base_pair, i, j));
}

template < class ValueT >
inline
ValueT&
ParameterHash<ValueT>::
base_pair(NUCL i, NUCL j)
{
  return get_by_key(string_format(format_base_pair, i, j));
}
#endif

#if PARAMS_BASE_PAIR_DIST
static const char* format_base_pair_dist_at_least = "base_pair_dist_at_least_%d";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
base_pair_dist_at_least(uint l) const
{
  return get_by_key(string_format(format_base_pair_dist_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
base_pair_dist_at_least(uint l)
{
  return get_by_key(string_format(format_base_pair_dist_at_least, l));
}
#endif

#if PARAMS_TERMINAL_MISMATCH
static const char* format_terminal_mismatch = "terminal_mismatch_%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
  return get_by_key(string_format(format_terminal_mismatch, i1, j1, i2, j2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2)
{
  return get_by_key(string_format(format_terminal_mismatch, i1, j1, i2, j2));
}
#endif

#if PARAMS_HAIRPIN_LENGTH
static const char* format_hairpin_length_at_least = "hairpin_length_at_least_%d";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_length_at_least(uint l) const
{
  return get_by_key(string_format(format_hairpin_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_length_at_least(uint l)
{
  return get_by_key(string_format(format_hairpin_length_at_least, l));
}
#endif

static const std::string format_hairpin_nucleotides("hairpin_nucleotides_");

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l) const
{
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
  std::string h(l, ' ');
  auto x = std::copy(&s[i], &s[i+l], h.begin());
  return get_by_key(format_hairpin_nucleotides + h);
}

template <class ValueT>
inline
std::vector<int>
ParameterHash<ValueT>::
hairpin_nucleotides_cache(const std::vector<NUCL>& s, uint i, uint max_l) const
{
  std::vector<int> v(max_l, trie_t::CEDAR_NO_VALUE);
  size_t node_pos=0, key_pos=0;
  auto r = trie_.traverse(format_hairpin_nucleotides.c_str(), node_pos, key_pos, 
                          format_hairpin_nucleotides.size());
  if (r == trie_t::CEDAR_NO_PATH) return v;

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
hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l, const std::vector<int>& pos) const
{
  return pos[l]>=0 ? values_[pos[l]] : static_cast<ValueT>(0);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l, const std::vector<int>& pos)
{
  return pos[l]>=0 ? values_[pos[l]] : hairpin_nucleotides(s, i, l);
}

#if PARAMS_HELIX_LENGTH
static const char* format_helix_length_at_least = "helix_length_at_least_%d";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
helix_length_at_least(uint l) const
{
  return get_by_key(string_format(format_helix_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
helix_length_at_least(uint l)
{
  return get_by_key(string_format(format_helix_length_at_least, l));
}
#endif

#if PARAMS_ISOLATED_BASE_PAIR
static const char* format_isolated_base_pair = "isolated_base_pair";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
isolated_base_pair() const
{
  return get_by_key(format_isolated_base_pair);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
isolated_base_pair()
{
  return get_by_key(format_isolated_base_pair);
}
#endif

#if PARAMS_INTERNAL_EXPLICIT
static const char* format_internal_explicit = "internal_explicit_%d_%d";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_explicit(uint i, uint j) const
{
  return get_by_key(string_format(format_internal_explicit, i, j));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_explicit(uint i, uint j)
{
  return get_by_key(string_format(format_internal_explicit, i, j));
}
#endif

#if PARAMS_BULGE_LENGTH
static const char* format_bulge_length_at_least = "bulge_length_at_least_%d";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_length_at_least(uint l) const
{
  return get_by_key(string_format(format_bulge_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_length_at_least(uint l)
{
  return get_by_key(string_format(format_bulge_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_LENGTH
static const char* format_internal_length_at_least = "internal_length_at_least_%d";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_length_at_least(uint l) const
{
  return get_by_key(string_format(format_internal_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_length_at_least(uint l)
{
  return get_by_key(string_format(format_internal_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_SYMMETRY
static const char* format_internal_symmetric_length_at_least = "internal_symmetric_length_at_least_%d";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_symmetric_length_at_least(uint l) const
{
  return get_by_key(string_format(format_internal_symmetric_length_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_symmetric_length_at_least(uint l)
{
  return get_by_key(string_format(format_internal_symmetric_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_ASYMMETRY
static const char* format_internal_asymmetry_at_least = "internal_asymmetry_at_least_%d";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_asymmetry_at_least(uint l) const
{
  return get_by_key(string_format(format_internal_asymmetry_at_least, l));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_asymmetry_at_least(uint l)
{
  return get_by_key(string_format(format_internal_asymmetry_at_least, l));
}
#endif

static const std::string format_internal_nucleotides("internal_nucleotides_");

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m) const
{
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
  std::string nuc1(l+m+1, ' ');
  auto x = std::copy(&s[i], &s[i+l], nuc1.begin());
  *(x++) = '_';
  std::reverse_copy(&s[j-m+1], &s[j+1], x);

  auto key1 = format_internal_nucleotides+nuc1;
  auto k = trie_.template exactMatchSearch<int>(key1.c_str());
  if (k!=trie_t::CEDAR_NO_VALUE)
    return values_[k];

  trie_.update(key1.c_str()) = values_.size();

  std::string nuc2(l+m+1, ' ');
  x = std::copy(&s[j-m+1], &s[j+1], nuc2.begin());
  *(x++) = '_';
  std::reverse_copy(&s[i], &s[i+l], x);
  auto key2 = format_internal_nucleotides+nuc2;
  trie_.update(key2.c_str()) = values_.size();

  values_.push_back(static_cast<ValueT>(0.0));
  return values_.back();
}

template <class ValueT>
std::vector<std::vector<int>>
ParameterHash<ValueT>::
internal_nucleotides_cache(const std::vector<NUCL>& s, uint i, uint j,
                           uint max_l, uint max_m) const
{
  std::vector<std::vector<int>> ret(max_l+1, std::vector<int>(max_m+1, trie_t::CEDAR_NO_VALUE));

  size_t node_pos=0, key_pos=0;
  auto k = trie_.traverse(format_internal_nucleotides.c_str(), node_pos, key_pos, 
                          format_internal_nucleotides.size());

  for (uint l=0; l<=max_l; ++l)
  {
    if (l>0)
    {
      key_pos = 0;
      k = trie_.traverse(&s[i+l-1], node_pos, key_pos, 1);
      if (k==trie_t::CEDAR_NO_PATH) break;
    }

    auto node_pos2 = node_pos;
    key_pos = 0;
    k = trie_.traverse("_", node_pos2, key_pos, 1);

    for (uint m=0; m<=max_m && l+m<=DEFAULT_C_MAX_SINGLE_LENGTH; ++m)
    {
      if (l+m<1) continue;
      if (m>0)
      {
        key_pos = 0;
        k = trie_.traverse(&s[j-m+1], node_pos2, key_pos, 1);
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
internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m,
                     const std::vector<std::vector<int>>& pos) const
{
#if 0
#ifndef NDEBUG
  auto v = pos[l][m]>=0 ? values_[pos[l][m]] : static_cast<ValueT>(0);
  //std::cout << s[i] << " " << i << " " << l << " " << j << " " << m << " " << pos[l][m] << " " << v << " " << internal_nucleotides(s, i, l, j, m)<< std::endl;
  assert(v == internal_nucleotides(s, i, l, j, m));
#endif
#endif
  return pos[l][m]>=0 ? values_[pos[l][m]] : static_cast<ValueT>(0);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m,
                     const std::vector<std::vector<int>>& pos)
{
  return pos[l][m]>=0 ? values_[pos[l][m]] : internal_nucleotides(s, i, l, j, m);
}



#if PARAMS_HELIX_STACKING
static const char* format_helix_stacking = "helix_stacking_%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
  return get_by_key(string_format(format_helix_stacking, i1, j1, i2, j2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2)
{
  return get_by_key(string_format(format_helix_stacking, i1, j1, i2, j2));
}
#endif

#if PARAMS_HELIX_CLOSING
static const char* format_helix_closing = "helix_closing_%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
helix_closing(NUCL i, NUCL j) const
{
  return get_by_key(string_format(format_helix_closing, i, j));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
helix_closing(NUCL i, NUCL j)
{
  return get_by_key(string_format(format_helix_closing, i, j));
}
#endif

#if PARAMS_MULTI_LENGTH
static const char* format_multi_base = "multi_base";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
multi_base() const
{
  return get_by_key(format_multi_base);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
multi_base()
{
  return get_by_key(format_multi_base);
}

static const char* format_multi_unpaired = "multi_unpaired";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
multi_unpaired() const
{
  return get_by_key(format_multi_unpaired);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
multi_unpaired()
{
  return get_by_key(format_multi_unpaired);
}

static const char* format_multi_paired = "multi_paired";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
multi_paired() const
{
  return get_by_key(format_multi_paired);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
multi_paired()
{
  return get_by_key(format_multi_paired);
}
#endif

#if PARAMS_DANGLE
static const char* format_dangle_left = "dangle_left_%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
dangle_left(NUCL i1, NUCL j1, NUCL i2) const
{
  return get_by_key(string_format(format_dangle_left, i1, j1, i2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
dangle_left(NUCL i1, NUCL j1, NUCL i2)
{
  return get_by_key(string_format(format_dangle_left, i1, j1, i2));
}

static const char* format_dangle_right = "dangle_right_%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
dangle_right(NUCL i1, NUCL j1, NUCL j2) const
{
  return get_by_key(string_format(format_dangle_right, i1, j1, j2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
dangle_right(NUCL i1, NUCL j1, NUCL j2)
{
  return get_by_key(string_format(format_dangle_right, i1, j1, j2));
}
#endif

#if PARAMS_EXTERNAL_LENGTH
static const char* format_external_unpaired = "external_unpaired";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
external_unpaired() const
{
  return get_by_key(format_external_unpaired);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
external_unpaired()
{
  return get_by_key(format_external_unpaired);
}

static const char* format_external_paired = "external_paired";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
external_paired() const
{
  return get_by_key(format_external_paired);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
external_paired()
{
  return get_by_key(format_external_paired);
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
