#include "../config.h"
#include "FeatureMap.hpp"

FeatureMap::
FeatureMap(const char* def_bases /*= "ACGUP" */)
  : def_bases_(def_bases), NBASES(def_bases_.size()), hash_(), keys_()
#ifdef USE_CACHE
#if PARAMS_BASE_PAIR
  , cache_base_pair_(NBASES, VI(NBASES, -1))
#endif
#if PARAMS_BASE_PAIR_DIST
  , cache_base_pair_dist_at_least_(D_MAX_BP_DIST_THRESHOLDS, -1)
#endif
#if PARAMS_TERMINAL_MISMATCH
  , cache_terminal_mismatch_(NBASES, VVVI(NBASES, VVI(NBASES, VI(NBASES, -1))))
#endif
#if PARAMS_HAIRPIN_LENGTH
  , cache_hairpin_length_at_least_(D_MAX_HAIRPIN_LENGTH, -1)
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
#if PARAMS_HELIX_STACKING
  , cache_helix_stacking_(NBASES, VVVI(NBASES, VVI(NBASES, VI(NBASES, -1))))
#endif
#if PARAMS_HELIX_CLOSING
  , cache_helix_closing_(NBASES, VI(NBASES, -1))
#endif
#if PARAMS_MULTI_LENGTH
  , cache_multi_base_(-1)
  , cache_multi_unpaired_(-1)
  , cache_multi_paired_(-1)
#endif
#if PARAMS_DANGLE
  , cache_dangle_left_(NBASES, VVI(NBASES, VI(NBASES, -1)))
  , cache_dangle_right_(NBASES, VVI(NBASES, VI(NBASES, -1)))
#endif
#if PARAMS_EXTERNAL_LENGTH
  , cache_external_unpaired_(-1)
  , cache_external_paired_(-1)
#endif
#endif
{
  std::fill(std::begin(is_base_), std::end(is_base_), -1);
  for (size_t i=0; i!=def_bases_.size(); ++i)
    is_base_[def_bases[i]] = i;

  initialize_cache();
}

void
FeatureMap::
initialize_cache()
{
#ifdef USE_CACHE
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
#endif
}

inline
size_t
FeatureMap::
set_key(const std::string& key)
{
  auto r = hash_.emplace(key, keys_.size());
  if (!r.second) return r.first->second;
  auto k = key;

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
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }
  else if (key.find(s_base_pair) == 0 && key.size() == s_base_pair.size()+2)
  {
    std::string nuc1 = key.substr(s_base_pair.size()+0, 1);
    std::string nuc2 = key.substr(s_base_pair.size()+1, 1);
    std::string key2 = s_base_pair + nuc2 + nuc1;
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }
  else if (key.find(s_helix_stacking) == 0)
  {
    std::string nuc1 = key.substr(s_helix_stacking.size()+0, 1);
    std::string nuc2 = key.substr(s_helix_stacking.size()+1, 1);
    std::string nuc3 = key.substr(s_helix_stacking.size()+2, 1);
    std::string nuc4 = key.substr(s_helix_stacking.size()+3, 1);
    std::string key2 = s_helix_stacking + nuc4 + nuc3 + nuc2 + nuc1;
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }
  else if (key.find(s_internal_explicit) == 0)
  {
    size_t pos=key.find("_", s_internal_explicit.size());
    std::string l1 = key.substr(s_internal_explicit.size(), pos-s_internal_explicit.size());
    std::string l2 = key.substr(pos+1);
    std::string key2 = s_internal_explicit + l2 + "_" + l1;
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }

  keys_.push_back(k.c_str());
  return keys_.size()-1;
}


#if PARAMS_BASE_PAIR
static const char* format_base_pair = "base_pair_%c%c";

void
FeatureMap::
initialize_cache_base_pair()
{
#ifdef USE_CACHE
  for (auto i : def_bases_)
    for (auto j : def_bases_)
      cache_base_pair_[is_base_[i]][is_base_[j]] = set_key(SPrintF(format_base_pair, i, j));
#endif
}

inline
size_t
FeatureMap::
base_pair(NUCL i, NUCL j) const
{
#ifdef USE_CACHE
  auto ii=is_base_[i], jj=is_base_[j];
  if (ii>=0 && jj>=0)
    return cache_base_pair_[ii][jj];
#endif
  return find_key(std::string("base_pair_")+i+j);
}
#endif

#if PARAMS_BASE_PAIR_DIST
static const char* format_base_pair_dist_at_least = "base_pair_dist_at_least_%d";

void
FeatureMap::
initialize_cache_base_pair_dist_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_base_pair_dist_at_least_.size(); ++i)
    cache_base_pair_dist_at_least_[i] = set_key(SPrintF(format_base_pair_dist_at_least, i));
#endif
}

inline
size_t
FeatureMap::
base_pair_dist_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_base_pair_dist_at_least_.size())
    return cache_base_pair_dist_at_least_[l];
#endif
  return find_key(SPrintF(format_base_pair_dist_at_least, l));
}
#endif

#if PARAMS_TERMINAL_MISMATCH
static const char* format_terminal_mismatch = "terminal_mismatch_%c%c%c%c";

void
FeatureMap::
initialize_cache_terminal_mismatch()
{
#ifdef USE_CACHE
  for (auto i1: def_bases_)
  {
    auto ii1 = is_base_[i1];
    for (auto j1: def_bases_)
    {
      auto jj1 = is_base_[j1];
      for (auto i2: def_bases_)
      {
        auto ii2 = is_base_[i2];
        for (auto j2: def_bases_)
        {
          auto jj2 = is_base_[j2];
          cache_terminal_mismatch_[ii1][jj1][ii2][jj2]
            = set_key(SPrintF(format_terminal_mismatch, i1, j1, i2, j2));
        }
      }
    }
  }
#endif
}

inline
size_t
FeatureMap::
terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
#ifdef USE_CACHE
  auto ii1 = is_base_[i1], jj1 = is_base_[j1];
  auto ii2 = is_base_[i2], jj2 = is_base_[j2];
  if (ii1>=0 && jj1>=0 && ii2>=0 && jj2>=0)
    return cache_terminal_mismatch_[ii1][jj1][ii2][jj2];
#endif
  return find_key(SPrintF(format_terminal_mismatch, i1, j1, i2, j2));
}
#endif

#if PARAMS_HAIRPIN_LENGTH
static const char* format_hairpin_length_at_least = "hairpin_length_at_least_%d";

void
FeatureMap::
initialize_cache_hairpin_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_hairpin_length_at_least_.size(); ++i)
    cache_hairpin_length_at_least_[i] = set_key(SPrintF(format_hairpin_length_at_least, i));
#endif
}

inline
size_t
FeatureMap::
hairpin_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_hairpin_length_at_least_.size())
    return cache_hairpin_length_at_least_[l];
#endif
  return find_key(SPrintF(format_hairpin_length_at_least, l));
}
#endif

static const std::string format_hairpin_nucleotides("hairpin_nucleotides_");

inline
size_t
FeatureMap::
hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l) const
{
  std::string h(l, ' ');
  std::copy(&s[i], &s[i]+l, h.begin());
  return find_key(format_hairpin_nucleotides + h);
}

#if PARAMS_HELIX_LENGTH
static const char* format_helix_length_at_least = "helix_length_at_least_%d";

void
FeatureMap::
initialize_cache_helix_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_helix_length_at_least_.size(); ++i)
    cache_helix_length_at_least_[i] = set_key(SPrintF(format_helix_length_at_least, i));
#endif
}

inline
size_t
FeatureMap::
helix_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_helix_length_at_least_.size())
    return cache_helix_length_at_least_[l];
#endif
  return find_key(SPrintF(format_helix_length_at_least, l));
}
#endif

#if PARAMS_ISOLATED_BASE_PAIR
static const char* format_isolated_base_pair = "isolated_base_pair";

void
FeatureMap::
initialize_cache_isolated_base_pair()
{
#ifdef USE_CACHE
  cache_isolated_base_pair_ = set_key(format_isolated_base_pair);
#endif
}

inline
size_t
FeatureMap::
isolated_base_pair() const
{
#ifdef USE_CACHE
  return cache_isolated_base_pair_;
#else
  return find_key(format_isolated_base_pair);
#endif
}
#endif

#if PARAMS_INTERNAL_EXPLICIT
static const char* format_internal_explicit = "internal_explicit_%d_%d";

void
FeatureMap::
initialize_cache_internal_explicit()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_internal_explicit_.size(); ++i)
    for (size_t j=0; j!=cache_internal_explicit_[i].size(); ++j)
      cache_internal_explicit_[i][j] = set_key(SPrintF(format_internal_explicit, i, j));
#endif
}

inline
size_t
FeatureMap::
internal_explicit(uint i, uint j) const
{
#ifdef USE_CACHE
  if (i<cache_internal_explicit_.size() && j<cache_internal_explicit_[i].size())
    return cache_internal_explicit_[i][j];
#endif
  return find_key(SPrintF(format_internal_explicit, i, j));
}
#endif

#if PARAMS_BULGE_LENGTH
static const char* format_bulge_length_at_least = "bulge_length_at_least_%d";

void
FeatureMap::
initialize_cache_bulge_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_bulge_length_at_least_.size(); ++i)
    cache_bulge_length_at_least_[i] = set_key(SPrintF(format_bulge_length_at_least, i));
#endif
}

inline
size_t
FeatureMap::
bulge_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_bulge_length_at_least_.size())
    return cache_bulge_length_at_least_[l];
#endif
  return find_key(SPrintF(format_bulge_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_LENGTH
static const char* format_internal_length_at_least = "internal_length_at_least_%d";

void
FeatureMap::
initialize_cache_internal_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_internal_length_at_least_.size(); ++i)
    cache_internal_length_at_least_[i] = set_key(SPrintF(format_internal_length_at_least, i));
#endif
}

inline
size_t
FeatureMap::
internal_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_internal_length_at_least_.size())
    return cache_internal_length_at_least_[l];
#endif
  return find_key(SPrintF(format_internal_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_SYMMETRY
static const char* format_internal_symmetric_length_at_least = "internal_symmetric_length_at_least_%d";

void
FeatureMap::
initialize_cache_internal_symmetric_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_internal_symmetric_length_at_least_.size(); ++i)
    cache_internal_symmetric_length_at_least_[i] = set_key(SPrintF(format_internal_symmetric_length_at_least, i));
#endif
}

inline
size_t
FeatureMap::
internal_symmetric_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_internal_symmetric_length_at_least_.size())
    return cache_internal_symmetric_length_at_least_[l];
#endif
  return find_key(SPrintF(format_internal_symmetric_length_at_least, l));
}
#endif

#if PARAMS_INTERNAL_ASYMMETRY
static const char* format_internal_asymmetry_at_least = "internal_asymmetry_at_least_%d";

void
FeatureMap::
initialize_cache_internal_asymmetry_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_internal_asymmetry_at_least_.size(); ++i)
    cache_internal_asymmetry_at_least_[i] = set_key(SPrintF(format_internal_asymmetry_at_least, i));
#endif
}

inline
size_t
FeatureMap::
internal_asymmetry_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_internal_asymmetry_at_least_.size())
    return cache_internal_asymmetry_at_least_[l];
#endif
  return find_key(SPrintF(format_internal_asymmetry_at_least, l));
}
#endif

static const std::string format_internal_nucleotides("internal_nucleotides_");

inline
size_t
FeatureMap::
internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m) const
{
  std::string nuc(l+m+1, ' ');
  auto x = std::copy(&s[i], &s[i+l], nuc.begin());
  *(x++) = '_';
  std::reverse_copy(&s[j-m+1], &s[j+1], x);
  return find_key(format_internal_nucleotides + nuc);
}

#if PARAMS_HELIX_STACKING
static const char* format_helix_stacking = "helix_stacking_%c%c%c%c";

void
FeatureMap::
initialize_cache_helix_stacking()
{
#ifdef USE_CACHE
  for (auto i1: def_bases_)
  {
    auto ii1 = is_base_[i1];
    for (auto j1: def_bases_)
    {
      auto jj1 = is_base_[j1];
      for (auto i2: def_bases_)
      {
        auto ii2 = is_base_[i2];
        for (auto j2: def_bases_)
        {
          auto jj2 = is_base_[j2];
          cache_helix_stacking_[ii1][jj1][ii2][jj2] = set_key(SPrintF(format_helix_stacking, i1, j1, i2, j2));
        }
      }
    }
  }
#endif
}

inline
size_t
FeatureMap::
helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
#ifdef USE_CACHE
  auto ii1=is_base_[i1], jj1=is_base_[j1];
  auto ii2=is_base_[i2], jj2=is_base_[j2];
  if (ii1>=0 && jj1>=0 && ii2>=0 && jj2>=0)
    return cache_helix_stacking_[ii1][jj1][ii2][jj2];
#endif
  return find_key(SPrintF(format_helix_stacking, i1, j1, i2, j2));
}
#endif

#if PARAMS_HELIX_CLOSING
static const char* format_helix_closing = "helix_closing_%c%c";

void
FeatureMap::
initialize_cache_helix_closing()
{
#ifdef USE_CACHE
  for (auto i: def_bases_)
  {
    auto ii = is_base_[i];
    for (auto j: def_bases_)
    {
      auto jj = is_base_[j];
      cache_helix_closing_[ii][jj] = set_key(SPrintF(format_helix_closing, i, j));
    }
  }
#endif
}

inline
size_t
FeatureMap::
helix_closing(NUCL i, NUCL j) const
{
#ifdef USE_CACHE
  auto ii=is_base_[i], jj=is_base_[j];
  if (ii>=0 && jj>=0)
    return cache_helix_closing_[ii][jj];
#endif
  return find_key(SPrintF(format_helix_closing, i, j));
}
#endif

#if PARAMS_MULTI_LENGTH
static const char* format_multi_base = "multi_base";

void
FeatureMap::
initialize_cache_multi_base()
{
#ifdef USE_CACHE
  cache_multi_base_ = set_key(format_multi_base);
#endif
}

inline
size_t
FeatureMap::
multi_base() const
{
#ifdef USE_CACHE
  return cache_multi_base_;
#else
  return find_key(format_multi_base);
#endif
}

static const char* format_multi_unpaired = "multi_unpaired";

void
FeatureMap::
initialize_cache_multi_unpaired()
{
#ifdef USE_CACHE
  cache_multi_unpaired_ = set_key(format_multi_unpaired);
#endif
}

inline
size_t
FeatureMap::
multi_unpaired() const
{
#ifdef USE_CACHE
  return cache_multi_unpaired_;
#else
  return find_key(format_multi_unpaired);
#endif
}

static const char* format_multi_paired = "multi_paired";

void
FeatureMap::
initialize_cache_multi_paired()
{
#ifdef USE_CACHE
  cache_multi_paired_ = set_key(format_multi_paired);
#endif
}

inline
size_t
FeatureMap::
multi_paired() const
{
#ifdef USE_CACHE
  return cache_multi_paired_;
#else
  return find_key(format_multi_paired);
#endif
}
#endif

#if PARAMS_DANGLE
static const char* format_dangle_left = "dangle_left_%c%c%c";

void
FeatureMap::
initialize_cache_dangle_left()
{
#ifdef USE_CACHE
  for (auto i1: def_bases_)
  {
    auto ii1 = is_base_[i1];
    for (auto j1: def_bases_)
    {
      auto jj1 = is_base_[j1];
      for (auto i2: def_bases_)
      {
        auto ii2 = is_base_[i2];
        cache_dangle_left_[ii1][jj1][ii2] = set_key(SPrintF(format_dangle_left, i1, j1, i2));
      }
    }
  }
#endif
}

inline
size_t
FeatureMap::
dangle_left(NUCL i1, NUCL j1, NUCL i2) const
{
#ifdef USE_CACHE
  auto ii1=is_base_[i1], jj1=is_base_[j1], ii2=is_base_[i2];
  if (ii1>=0 && jj1>=0 && ii2>=0)
    return cache_dangle_left_[ii1][jj1][ii2];
#endif
  return find_key(SPrintF(format_dangle_left, i1, j1, i2));
}

static const char* format_dangle_right = "dangle_right_%c%c%c";

inline
void
FeatureMap::
initialize_cache_dangle_right()
{
#ifdef USE_CACHE
  for (auto i1: def_bases_)
  {
    auto ii1 = is_base_[i1];
    for (auto j1: def_bases_)
    {
      auto jj1 = is_base_[j1];
      for (auto j2: def_bases_)
      {
        auto jj2 = is_base_[j2];
        cache_dangle_right_[ii1][jj1][jj2] = set_key(SPrintF(format_dangle_right, i1, j1, j2));
      }
    }
  }
#endif
}

inline
size_t
FeatureMap::
dangle_right(NUCL i1, NUCL j1, NUCL j2) const
{
#ifdef USE_CACHE
  auto ii1=is_base_[i1], jj1=is_base_[j1], jj2=is_base_[j2];
  if (ii1>=0 && jj1>=0 && jj2>=0)
    return cache_dangle_right_[ii1][jj1][jj2];
#endif
    return find_key(SPrintF(format_dangle_right, i1, j1, j2));
}
#endif

#if PARAMS_EXTERNAL_LENGTH
static const char* format_external_unpaired = "external_unpaired";

void
FeatureMap::
initialize_cache_external_unpaired()
{
#ifdef USE_CACHE
  cache_external_unpaired_ = set_key(format_external_unpaired);
#endif
}

inline
size_t
FeatureMap::
external_unpaired() const
{
#ifdef USE_CACHE
  return cache_external_unpaired_;
#else
  return find_key(format_external_unpaired);
#endif
}

static const char* format_external_paired = "external_paired";

void
FeatureMap::
initialize_cache_external_paired()
{
#ifdef USE_CACHE
  cache_external_paired_ = set_key(format_external_paired);
#endif
}

inline
size_t
FeatureMap::
external_paired() const
{
#ifdef USE_CACHE
  return cache_external_paired_;
#else
  return find_key(format_external_paired);
#endif
}
#endif
