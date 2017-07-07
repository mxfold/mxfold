#include "../config.h"
#include "Config.hpp"
#include "FeatureMap.hpp"
#include "LogSpace.hpp"
#include <sstream>

// constants
const std::string s_base_pair("base_pair_");
const std::string s_base_pair_dist_at_least("base_pair_dist_at_least_%d");
const std::string s_terminal_mismatch("terminal_mismatch_");
const std::string s_hairpin_length_at_least("hairpin_length_at_least_%d");
const std::string s_hairpin_nucleotides("hairpin_nucleotides_");
const std::string s_helix_length_at_least("helix_length_at_least_%d");
const std::string s_isolated_base_pair("isolated_base_pair");
const std::string s_internal_explicit("internal_explicit_%d_%d");
const std::string s_bulge_length_at_least("bulge_length_at_least_%d");
const std::string s_internal_length_at_least("internal_length_at_least_%d");
const std::string s_internal_symmetric_length_at_least("internal_symmetric_length_at_least_%d");
const std::string s_internal_asymmetry_at_least("internal_asymmetry_at_least_%d");
const std::string s_internal_nucleotides("internal_nucleotides_");
const std::string s_helix_stacking("helix_stacking_");
const std::string s_helix_closing("helix_closing_");
const std::string s_multi_base("multi_base");
const std::string s_multi_unpaired("multi_unpaired");
const std::string s_multi_paired("multi_paired");
const std::string s_dangle_left = "dangle_left_";
const std::string s_dangle_right = "dangle_right_";
const std::string s_external_unpaired = "external_unpaired";
const std::string s_external_paired = "external_paired";

// for vienna compatibility
const std::string s_internal_nucleotides_int11("internal_loop_11_");
const std::string s_internal_nucleotides_int21("internal_loop_21_");
const std::string s_internal_nucleotides_int12("internal_loop_12_");
const std::string s_internal_nucleotides_int22("internal_loop_22_");
const std::string s_terminal_mismatch_hairpin("terminal_mismatch_hairpin_");
const std::string s_terminal_mismatch_internal("terminal_mismatch_internal_");
const std::string s_terminal_mismatch_internal_1n("terminal_mismatch_internal_1n_");
const std::string s_terminal_mismatch_internal_23("terminal_mismatch_internal_23_");
const std::string s_terminal_mismatch_multi("terminal_mismatch_multi_");
const std::string s_terminal_mismatch_external("terminal_mismatch_external_");
const std::string s_ninio("ninio");
const std::string s_ninio_max("ninio_max");
const std::string s_termialAU("termialAU");
const std::string s_triloop("triloop_");
const std::string s_tetraloop("tetraloop_");
const std::string s_hexaloop("hexaloop_");

FeatureMap::
FeatureMap(const char* def_bases,
           const std::vector<std::string>& def_bps)
  : def_bases_(def_bases), NBASES(def_bases_.size()),
    def_bps_(def_bps), NBPS(def_bps_.size()),
    hash_(), keys_()
#ifdef USE_CACHE
#if PARAMS_BASE_PAIR
  , cache_base_pair_(NBPS, -1)
#endif
#if PARAMS_BASE_PAIR_DIST
  , cache_base_pair_dist_at_least_(D_MAX_BP_DIST_THRESHOLDS, -1)
#endif
#if PARAMS_TERMINAL_MISMATCH
  , cache_terminal_mismatch_(NBPS, VVI(NBASES, VI(NBASES, -1)))
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
  , cache_helix_stacking_(NBPS, VI(NBPS, -1))
#endif
#if PARAMS_HELIX_CLOSING
  , cache_helix_closing_(NBPS, -1)
#endif
#if PARAMS_MULTI_LENGTH
  , cache_multi_base_(-1)
  , cache_multi_unpaired_(-1)
  , cache_multi_paired_(-1)
#endif
#if PARAMS_DANGLE
  , cache_dangle_left_(NBPS, VI(NBASES, -1))
  , cache_dangle_right_(NBPS, VI(NBASES, -1))
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

  for (auto& e : is_complementary_)
    std::fill(std::begin(e), std::end(e), -1);
  for (size_t i=0; i!=def_bps_.size(); ++i)
    is_complementary_[def_bps_[i][0]][def_bps_[i][1]] = i;

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


size_t
FeatureMap::
find_key(const std::string& key) const
{
  auto itr = hash_.find(key);
  return itr != hash_.end() ? itr->second : -1u;
}

size_t
FeatureMap::
insert_key(const std::string& key)
{
  auto r = hash_.emplace(key, keys_.size());
  if (!r.second) return r.first->second;
  auto k = key;

  // symmetric features
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
    const char nuc1 = key[s_base_pair.size()+0];
    const char nuc2 = key[s_base_pair.size()+1];
    std::string key2 = s_base_pair + nuc2 + nuc1;
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }
  else if (key.find(s_helix_stacking) == 0)
  {
    const char nuc1 = key[s_helix_stacking.size()+0];
    const char nuc2 = key[s_helix_stacking.size()+1];
    const char nuc3 = key[s_helix_stacking.size()+2];
    const char nuc4 = key[s_helix_stacking.size()+3];
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
  else if (key.find(s_internal_nucleotides_int11) == 0)
  {
    const char nuc1 = key[s_internal_nucleotides_int11.size()+0];
    const char nuc2 = key[s_internal_nucleotides_int11.size()+1];
    const char nuc3 = key[s_internal_nucleotides_int11.size()+2];
    const char nuc4 = key[s_internal_nucleotides_int11.size()+3];
    const char nuc5 = key[s_internal_nucleotides_int11.size()+5];
    const char nuc6 = key[s_internal_nucleotides_int11.size()+7];
    std::string key2 = s_internal_nucleotides_int11 + nuc4 + nuc3 + nuc2 + nuc1 + "_" + nuc6 + "_" + nuc5;
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }
  else if (key.find(s_internal_nucleotides_int21) == 0)
  {
    const char nuc1 = key[s_internal_nucleotides_int21.size()+0];
    const char nuc2 = key[s_internal_nucleotides_int21.size()+1];
    const char nuc3 = key[s_internal_nucleotides_int21.size()+2];
    const char nuc4 = key[s_internal_nucleotides_int21.size()+3];
    const char nuc5 = key[s_internal_nucleotides_int21.size()+5];
    const char nuc6 = key[s_internal_nucleotides_int21.size()+6];
    const char nuc7 = key[s_internal_nucleotides_int21.size()+8];
    std::string key2 = s_internal_nucleotides_int12 + nuc4 + nuc3 + nuc2 + nuc1 + "_" + nuc7 + "_" + nuc6 + nuc5;
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }
  else if (key.find(s_internal_nucleotides_int12) == 0)
  {
    const char nuc1 = key[s_internal_nucleotides_int12.size()+0];
    const char nuc2 = key[s_internal_nucleotides_int12.size()+1];
    const char nuc3 = key[s_internal_nucleotides_int12.size()+2];
    const char nuc4 = key[s_internal_nucleotides_int12.size()+3];
    const char nuc5 = key[s_internal_nucleotides_int12.size()+5];
    const char nuc6 = key[s_internal_nucleotides_int12.size()+7];
    const char nuc7 = key[s_internal_nucleotides_int12.size()+8];
    std::string key2 = s_internal_nucleotides_int21 + nuc4 + nuc3 + nuc2 + nuc1 + "_" + nuc7 + nuc6 + "_" + nuc5;
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }
  else if (key.find(s_internal_nucleotides_int22) == 0)
  {
    const char nuc1 = key[s_internal_nucleotides_int22.size()+0];
    const char nuc2 = key[s_internal_nucleotides_int22.size()+1];
    const char nuc3 = key[s_internal_nucleotides_int22.size()+2];
    const char nuc4 = key[s_internal_nucleotides_int22.size()+3];
    const char nuc5 = key[s_internal_nucleotides_int22.size()+5];
    const char nuc6 = key[s_internal_nucleotides_int22.size()+6];
    const char nuc7 = key[s_internal_nucleotides_int22.size()+8];
    const char nuc8 = key[s_internal_nucleotides_int22.size()+9];
    std::string key2 = s_internal_nucleotides_int22 + nuc4 + nuc3 + nuc2 + nuc1 + "_" + nuc8 + nuc7 + "_" + nuc6 + nuc5;
    hash_[key2] = keys_.size();
    if (k>key2) k = key2;
  }

  keys_.push_back(k.c_str());
  return keys_.size()-1;
}

size_t
FeatureMap::
insert_keyval(const std::string& key, std::vector<param_value_type>& vals, param_value_type v)
{
  size_t i = insert_key(key);
  if (i>=vals.size()) vals.resize(i+1, 0.0);
  vals[i] = v;
  return i;
}

std::vector<param_value_type>
FeatureMap::
load_from_hash(const std::unordered_map<std::string, param_value_type>& h)
{
  std::vector<param_value_type> vals;
  hash_.clear();
  keys_.clear();
  for (auto e: h)
    insert_keyval(e.first, vals, e.second);

  initialize_cache()  ;
  return std::move(vals);
}

std::vector<param_value_type>
FeatureMap::
read_from_file(const std::string& filename)
{
  std::vector<param_value_type> vals;
  hash_.clear();
  keys_.clear();
  std::ifstream is(filename.c_str());
  if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  std::string k;
  param_value_type v;
  if (!std::isdigit(is.peek()))
  {
    while (is >> k >> v)
      if (v!=0.0)
        insert_keyval(k, vals, v);
  }
  else
  {
    // for reading AdaGradRDA outputs
    float eta, lambda, eps, t, s1, s2;
    is >> eta >> lambda >> eps >> t;
    while (is >> k >> v >> s1 >> s2)
      if (v!=0.0)
        insert_keyval(k, vals, v);
  }

  initialize_cache();
  return std::move(vals);
}

std::string
ignore_comment(const std::string& l)
{
  size_t cp1 = l.find("/*");
  if (cp1 == std::string::npos) return l;
  size_t cp2 = l.find("*/");
  if (cp2 == std::string::npos)
    throw std::runtime_error("unclosed comment");
  return l.substr(0, cp1)+l.substr(cp2+2);
}

std::vector<param_value_type>
get_array1(const std::string& l)
{
  std::vector<param_value_type> ret;

  std::istringstream is(ignore_comment(l));
  std::string s;
  while (is >> s)
  {
    if (s=="DEF")
      ret.push_back(0.);
    else if (s=="INF")
      ret.push_back(NEG_INF);
    else if (s=="NST")
      ret.push_back(50/100.);
    else
      ret.push_back(-atoi(s.c_str())/100.);
  }
  return std::move(ret);
}

std::vector<param_value_type>
get_array(std::istream& is, size_t sz)
{
  std::vector<param_value_type> v;
  v.reserve(sz);
  std::string l;
  while (v.size() < sz && std::getline(is, l))
  {
    auto w = get_array1(l);
    v.insert(v.end(), w.begin(), w.end());
  }
  return std::move(v);
}

std::pair<std::string, param_value_type>
get_loops(std::istream& is)
{
  std::string line;
  if (std::getline(is, line))
  {
    char loop[100];
    int v, w;
    if (sscanf(line.c_str(), "%s %d %d", loop, &v, &w)==3)
      return std::make_pair(std::string(loop), -v/100.);
  }
  return std::pair<std::string, param_value_type>();
}

std::vector<param_value_type>
FeatureMap::
import_from_vienna_parameters(const std::string& filename)
{
  std::vector<param_value_type> vals;
  hash_.clear();
  keys_.clear();
  std::ifstream is(filename.c_str());
  if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  std::string line;
  std::getline(is, line);
  if (line.compare(0, 30, "## RNAfold parameter file v2.0") != 0)
    throw std::runtime_error(std::string("Invalid file format: ") + filename);

  while (std::getline(is, line))
  {
    size_t pos = line.find("# ");
    if (pos == 0) {
      std::string ident = line.substr(pos+2);
      if (ident == "stack")
      {
        auto v = get_array(is, (NBPS+1)*(NBPS+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBPS+1; ++j, ++p)
            if (i<NBPS && j<NBPS)
              insert_keyval(s_helix_stacking+def_bps_[i]+def_bps_[j], vals, v[p]);
      }
      else if (ident == "mismatch_hairpin")
      {
        auto v = get_array(is, (NBPS+1)*(NBASES+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBASES+1; ++j)
            for (size_t k=0; k!=NBASES+1; ++k, ++p)
              if (i<NBPS && j!=0 && k!=0)
                insert_keyval(s_terminal_mismatch_hairpin+def_bps_[i]+def_bases_[j-1]+def_bases_[k-1], vals, v[p]);
      }
      else if (ident == "mismatch_interior")
      {
        auto v = get_array(is, (NBPS+1)*(NBASES+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBASES+1; ++j)
            for (size_t k=0; k!=NBASES+1; ++k, ++p)
              if (i<NBPS && j!=0 && k!=0)
                insert_keyval(s_terminal_mismatch_internal+def_bps_[i]+def_bases_[j-1]+def_bases_[k-1], vals, v[p]);
      }
      else if (ident == "mismatch_interior_1n")
      {
        auto v = get_array(is, (NBPS+1)*(NBASES+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBASES+1; ++j)
            for (size_t k=0; k!=NBASES+1; ++k, ++p)
              if (i<NBPS && j!=0 && k!=0)
                insert_keyval(s_terminal_mismatch_internal_1n+def_bps_[i]+def_bases_[j-1]+def_bases_[k-1], vals, v[p]);
      }
      else if (ident == "mismatch_interior_23")
      {
        auto v = get_array(is, (NBPS+1)*(NBASES+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBASES+1; ++j)
            for (size_t k=0; k!=NBASES+1; ++k, ++p)
              if (i<NBPS && j!=0 && k!=0)
                insert_keyval(s_terminal_mismatch_internal_23+def_bps_[i]+def_bases_[j-1]+def_bases_[k-1], vals, v[p]);
      }
      else if (ident == "mismatch_multi")
      {
        auto v = get_array(is, (NBPS+1)*(NBASES+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBASES+1; ++j)
            for (size_t k=0; k!=NBASES+1; ++k, ++p)
              if (i<NBPS && j!=0 && k!=0)
                insert_keyval(s_terminal_mismatch_multi+def_bps_[i]+def_bases_[j-1]+def_bases_[k-1], vals, v[p]);
      }
      else if (ident == "mismatch_exterior")
      {
        auto v = get_array(is, (NBPS+1)*(NBASES+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBASES+1; ++j)
            for (size_t k=0; k!=NBASES+1; ++k, ++p)
              if (i<NBPS && j!=0 && k!=0)
                insert_keyval(s_terminal_mismatch_external+def_bps_[i]+def_bases_[j-1]+def_bases_[k-1], vals, v[p]);
      }
      else if (ident == "dangle5")
      {
        auto v = get_array(is, (NBPS+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBASES+1; ++j, ++p)
            if (i<NBPS && j!=0)
              insert_keyval(s_dangle_left+def_bps_[i]+def_bases_[j-1], vals, v[p]);
      }
      else if (ident == "dangle3")
      {
        auto v = get_array(is, (NBPS+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBASES+1; ++j, ++p)
            if (i<NBPS && j!=0)
              insert_keyval(s_dangle_right+def_bps_[i]+def_bases_[j-1], vals, v[p]);
      }
      else if (ident == "int11")
      {
        auto v = get_array(is, (NBPS+1)*(NBPS+1)*(NBASES+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBPS+1; ++j)
            for (size_t k=0; k!=NBASES+1; ++k)
              for (size_t l=0; l!=NBASES+1; ++l, ++p)
                if (i<NBPS && j<NBPS && k!=0 && l!=0)
                  insert_keyval(s_internal_nucleotides_int11+def_bps_[i]+def_bps_[j]+"_"+def_bases_[k-1]+"_"+def_bases_[l-1], vals, v[p]);
      }
      else if (ident == "int21")
      {
        auto v = get_array(is, (NBPS+1)*(NBPS+1)*(NBASES+1)*(NBASES+1)*(NBASES+1));
        for (size_t i=0, p=0; i!=NBPS+1; ++i)
          for (size_t j=0; j!=NBPS+1; ++j)
            for (size_t k=0; k!=NBASES+1; ++k)
              for (size_t l=0; l!=NBASES+1; ++l)
                for (size_t m=0; m!=NBASES+1; ++m, ++p)
                if (i<NBPS && j<NBPS && k!=0 && l!=0 && m!=0)
                  insert_keyval(s_internal_nucleotides_int21+def_bps_[i]+def_bps_[j]+"_"+def_bases_[k-1]+def_bases_[l-1]+"_"+def_bases_[m-1], vals, v[p]);
      }
      else if (ident == "int22")
      {
        auto v = get_array(is, NBPS*NBPS*NBASES*NBASES*NBASES*NBASES);
        for (size_t i=0, p=0; i!=NBPS; ++i)
          for (size_t j=0; j!=NBPS; ++j)
            for (size_t k=0; k!=NBASES; ++k)
              for (size_t l=0; l!=NBASES; ++l)
                for (size_t m=0; m!=NBASES; ++m)
                  for (size_t n=0; n!=NBASES; ++n, ++p)
                    insert_keyval(s_internal_nucleotides_int22+def_bps_[i]+def_bps_[j]+"_"+def_bases_[k]+def_bases_[l]+"_"+def_bases_[m]+def_bases_[n], vals, v[p]);
      }
      else if (ident == "hairpin")
      {
        auto v = get_array(is, 31);
        for (size_t i = v.size()-1; i > 0; --i)
          insert_keyval(SPrintF(s_hairpin_length_at_least.c_str(), i), vals, v[i]-v[i-1]);
      }
      else if (ident == "bulge")
      {
        auto v = get_array(is, 31);
        for (size_t i = 0; i < v.size(); i++)
        {
          param_value_type w=0.0;
          for (size_t j = i; j < v.size(); j++) w += v[j];
          insert_keyval(SPrintF(s_bulge_length_at_least.c_str(), i), vals, w);
        }
      }
      else if (ident == "interior")
      {
        auto v = get_array(is, 31);
        for (size_t i = 0; i < v.size(); i++)
        {
          param_value_type w=0.0;
          for (size_t j = i; j < v.size(); j++) w += v[j];
          insert_keyval(SPrintF(s_internal_length_at_least.c_str(), i), vals, w);
        }
      }
      else if (ident == "NINIO")
      {
        auto v = get_array(is, 3);
        insert_keyval(s_ninio, vals, v[0]);
        insert_keyval(s_ninio_max, vals, v[2]);
      }
      else if (ident == "ML_params")
      {
        auto v = get_array(is, 6);
        insert_keyval(s_multi_unpaired, vals, v[0]);
        insert_keyval(s_multi_base, vals, v[2]);
        insert_keyval(s_multi_paired, vals, v[4]);
      }
      else if (ident == "Misc")
      {
        auto v = get_array(is, 6);
        insert_keyval(s_termialAU, vals, v[2]);
      }
      else if (ident == "Triloops")
      {
        while (1)
        {
          auto v = get_loops(is);
          if (v.first.empty()) break;
          insert_keyval(s_triloop+v.first, vals, v.second);
        }
      }
      else if (ident == "Tetraloops")
      {
        while (1)
        {
          auto v = get_loops(is);
          if (v.first.empty()) break;
          insert_keyval(s_tetraloop+v.first, vals, v.second);
        }
      }
      else if (ident =="Hexaloops")
      {
        while (1)
        {
          auto v = get_loops(is);
          if (v.first.empty()) break;
          insert_keyval(s_hexaloop+v.first, vals, v.second);
        }
      }
    }
  }

  //initialize_cache();
  return std::move(vals);
}

void
FeatureMap::
write_to_file(const std::string& filename, const std::vector<param_value_type>& vals) const
{
  std::ofstream os(filename.c_str());
  if (!os) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  std::vector<size_t> idx(keys_.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&](size_t i, size_t j) { return keys_[i] < keys_[j]; });
  for (auto i: idx)
    if (vals[i]!=0.0)
      os << keys_[i] << " " << vals[i] << std::endl;
}

#if PARAMS_BASE_PAIR
void
FeatureMap::
initialize_cache_base_pair()
{
#ifdef USE_CACHE
  for (auto bp : def_bps_)
    cache_base_pair_[is_complementary_[bp[0]][bp[1]]] = insert_key(s_base_pair+bp[0]+bp[1]);
#endif
}

size_t
FeatureMap::
find_base_pair(NUCL i, NUCL j) const
{
#ifdef USE_CACHE
  auto ij=is_complementary_[i][j];
  if (ij>=0)
    return cache_base_pair_[ij];
#endif
  return find_key(s_base_pair+i+j);
}

size_t
FeatureMap::
insert_base_pair(NUCL i, NUCL j)
{
#ifdef USE_CACHE
  auto ij=is_complementary_[i][j];
  if (ij>=0)
    return cache_base_pair_[ij];
#endif
  return insert_key(s_base_pair+i+j);
}
#endif

#if PARAMS_BASE_PAIR_DIST
void
FeatureMap::
initialize_cache_base_pair_dist_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_base_pair_dist_at_least_.size(); ++i)
    cache_base_pair_dist_at_least_[i] = insert_key(SPrintF(s_base_pair_dist_at_least.c_str(), i));
#endif
}

size_t
FeatureMap::
find_base_pair_dist_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_base_pair_dist_at_least_.size())
    return cache_base_pair_dist_at_least_[l];
#endif
  return find_key(SPrintF(s_base_pair_dist_at_least.c_str(), l));
}

size_t
FeatureMap::
insert_base_pair_dist_at_least(uint l)
{
#ifdef USE_CACHE
  if (l<cache_base_pair_dist_at_least_.size())
    return cache_base_pair_dist_at_least_[l];
#endif
  return insert_key(SPrintF(s_base_pair_dist_at_least.c_str(), l));
}
#endif

#if PARAMS_TERMINAL_MISMATCH
void
FeatureMap::
initialize_cache_terminal_mismatch()
{
#ifdef USE_CACHE
  for (auto bp : def_bps_)
  {
    auto ij1 = is_complementary_[bp[0]][bp[1]];
    for (auto i2: def_bases_)
    {
      auto ii2 = is_base_[i2];
      for (auto j2: def_bases_)
      {
        auto jj2 = is_base_[j2];
        cache_terminal_mismatch_[ij1][ii2][jj2] = insert_key(s_terminal_mismatch+bp+i2+j2);
      }
    }
  }
#endif
}

size_t
FeatureMap::
find_terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
#ifdef USE_CACHE
  auto ij1 = is_complementary_[i1][j1];
  auto ii2 = is_base_[i2], jj2 = is_base_[j2];
  if (ij1>=0 && ii2>=0 && jj2>=0)
    return cache_terminal_mismatch_[ij1][ii2][jj2];
#endif
  return find_key(s_terminal_mismatch+i1+j1+i2+j2);
}

size_t
FeatureMap::
insert_terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2)
{
#ifdef USE_CACHE
  auto ij1 = is_complementary_[i1][j1];
  auto ii2 = is_base_[i2], jj2 = is_base_[j2];
  if (ij1>=0 && ii2>=0 && jj2>=0)
    return cache_terminal_mismatch_[ij1][ii2][jj2];
#endif
  return insert_key(s_terminal_mismatch+i1+j1+i2+j2);
}
#endif

#if PARAMS_HAIRPIN_LENGTH
void
FeatureMap::
initialize_cache_hairpin_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_hairpin_length_at_least_.size(); ++i)
    cache_hairpin_length_at_least_[i] = insert_key(SPrintF(s_hairpin_length_at_least.c_str(), i));
#endif
}

size_t
FeatureMap::
find_hairpin_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_hairpin_length_at_least_.size())
    return cache_hairpin_length_at_least_[l];
#endif
  return find_key(SPrintF(s_hairpin_length_at_least.c_str(), l));
}

size_t
FeatureMap::
insert_hairpin_length_at_least(uint l)
{
#ifdef USE_CACHE
  if (l<cache_hairpin_length_at_least_.size())
    return cache_hairpin_length_at_least_[l];
#endif
  return insert_key(SPrintF(s_hairpin_length_at_least.c_str(), l));
}
#endif

size_t
FeatureMap::
find_hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l) const
{
  std::string h(l, ' ');
  std::copy(&s[i], &s[i]+l, h.begin());
  return find_key(s_hairpin_nucleotides + h);
}

size_t
FeatureMap::
insert_hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l)
{
  std::string h(l, ' ');
  std::copy(&s[i], &s[i]+l, h.begin());
  return insert_key(s_hairpin_nucleotides + h);
}

#if PARAMS_HELIX_LENGTH
void
FeatureMap::
initialize_cache_helix_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_helix_length_at_least_.size(); ++i)
    cache_helix_length_at_least_[i] = insert_key(SPrintF(s_helix_length_at_least.c_str(), i));
#endif
}

size_t
FeatureMap::
find_helix_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_helix_length_at_least_.size())
    return cache_helix_length_at_least_[l];
#endif
  return find_key(SPrintF(s_helix_length_at_least.c_str(), l));
}

size_t
FeatureMap::
insert_helix_length_at_least(uint l)
{
#ifdef USE_CACHE
  if (l<cache_helix_length_at_least_.size())
    return cache_helix_length_at_least_[l];
#endif
  return insert_key(SPrintF(s_helix_length_at_least.c_str(), l));
}
#endif

#if PARAMS_ISOLATED_BASE_PAIR
void
FeatureMap::
initialize_cache_isolated_base_pair()
{
#ifdef USE_CACHE
  cache_isolated_base_pair_ = insert_key(s_isolated_base_pair);
#endif
}

size_t
FeatureMap::
find_isolated_base_pair() const
{
#ifdef USE_CACHE
  return cache_isolated_base_pair_;
#else
  return find_key(s_isolated_base_pair);
#endif
}

size_t
FeatureMap::
insert_isolated_base_pair()
{
#ifdef USE_CACHE
  return cache_isolated_base_pair_;
#else
  return insert_key(s_isolated_base_pair);
#endif
}
#endif

#if PARAMS_INTERNAL_EXPLICIT
void
FeatureMap::
initialize_cache_internal_explicit()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_internal_explicit_.size(); ++i)
    for (size_t j=0; j!=cache_internal_explicit_[i].size(); ++j)
      cache_internal_explicit_[i][j] = insert_key(SPrintF(s_internal_explicit.c_str(), i, j));
#endif
}

size_t
FeatureMap::
find_internal_explicit(uint i, uint j) const
{
#ifdef USE_CACHE
  if (i<cache_internal_explicit_.size() && j<cache_internal_explicit_[i].size())
    return cache_internal_explicit_[i][j];
#endif
  return find_key(SPrintF(s_internal_explicit.c_str(), i, j));
}

size_t
FeatureMap::
insert_internal_explicit(uint i, uint j)
{
#ifdef USE_CACHE
  if (i<cache_internal_explicit_.size() && j<cache_internal_explicit_[i].size())
    return cache_internal_explicit_[i][j];
#endif
  return insert_key(SPrintF(s_internal_explicit.c_str(), i, j));
}
#endif

#if PARAMS_BULGE_LENGTH
void
FeatureMap::
initialize_cache_bulge_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_bulge_length_at_least_.size(); ++i)
    cache_bulge_length_at_least_[i] = insert_key(SPrintF(s_bulge_length_at_least.c_str(), i));
#endif
}

size_t
FeatureMap::
find_bulge_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_bulge_length_at_least_.size())
    return cache_bulge_length_at_least_[l];
#endif
  return find_key(SPrintF(s_bulge_length_at_least.c_str(), l));
}

size_t
FeatureMap::
insert_bulge_length_at_least(uint l)
{
#ifdef USE_CACHE
  if (l<cache_bulge_length_at_least_.size())
    return cache_bulge_length_at_least_[l];
#endif
  return insert_key(SPrintF(s_bulge_length_at_least.c_str(), l));
}
#endif

#if PARAMS_INTERNAL_LENGTH
void
FeatureMap::
initialize_cache_internal_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_internal_length_at_least_.size(); ++i)
    cache_internal_length_at_least_[i] = insert_key(SPrintF(s_internal_length_at_least.c_str(), i));
#endif
}

size_t
FeatureMap::
find_internal_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_internal_length_at_least_.size())
    return cache_internal_length_at_least_[l];
#endif
  return find_key(SPrintF(s_internal_length_at_least.c_str(), l));
}

size_t
FeatureMap::
insert_internal_length_at_least(uint l)
{
#ifdef USE_CACHE
  if (l<cache_internal_length_at_least_.size())
    return cache_internal_length_at_least_[l];
#endif
  return insert_key(SPrintF(s_internal_length_at_least.c_str(), l));
}
#endif

#if PARAMS_INTERNAL_SYMMETRY
void
FeatureMap::
initialize_cache_internal_symmetric_length_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_internal_symmetric_length_at_least_.size(); ++i)
    cache_internal_symmetric_length_at_least_[i] = insert_key(SPrintF(s_internal_symmetric_length_at_least.c_str(), i));
#endif
}

size_t
FeatureMap::
find_internal_symmetric_length_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_internal_symmetric_length_at_least_.size())
    return cache_internal_symmetric_length_at_least_[l];
#endif
  return find_key(SPrintF(s_internal_symmetric_length_at_least.c_str(), l));
}

size_t
FeatureMap::
insert_internal_symmetric_length_at_least(uint l)
{
#ifdef USE_CACHE
  if (l<cache_internal_symmetric_length_at_least_.size())
    return cache_internal_symmetric_length_at_least_[l];
#endif
  return insert_key(SPrintF(s_internal_symmetric_length_at_least.c_str(), l));
}
#endif

#if PARAMS_INTERNAL_ASYMMETRY
void
FeatureMap::
initialize_cache_internal_asymmetry_at_least()
{
#ifdef USE_CACHE
  for (size_t i=0; i!=cache_internal_asymmetry_at_least_.size(); ++i)
    cache_internal_asymmetry_at_least_[i] = insert_key(SPrintF(s_internal_asymmetry_at_least.c_str(), i));
#endif
}

size_t
FeatureMap::
find_internal_asymmetry_at_least(uint l) const
{
#ifdef USE_CACHE
  if (l<cache_internal_asymmetry_at_least_.size())
    return cache_internal_asymmetry_at_least_[l];
#endif
  return find_key(SPrintF(s_internal_asymmetry_at_least.c_str(), l));
}

size_t
FeatureMap::
insert_internal_asymmetry_at_least(uint l)
{
#ifdef USE_CACHE
  if (l<cache_internal_asymmetry_at_least_.size())
    return cache_internal_asymmetry_at_least_[l];
#endif
  return insert_key(SPrintF(s_internal_asymmetry_at_least.c_str(), l));
}
#endif

size_t
FeatureMap::
find_internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m) const
{
  std::string nuc(l+m+1, ' ');
  auto x = std::copy(&s[i], &s[i+l], nuc.begin());
  *(x++) = '_';
  std::reverse_copy(&s[j-m+1], &s[j+1], x);
  return find_key(s_internal_nucleotides + nuc);
}

size_t
FeatureMap::
insert_internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m)
{
  std::string nuc(l+m+1, ' ');
  auto x = std::copy(&s[i], &s[i+l], nuc.begin());
  *(x++) = '_';
  std::reverse_copy(&s[j-m+1], &s[j+1], x);
  return insert_key(s_internal_nucleotides + nuc);
}

#if PARAMS_HELIX_STACKING
void
FeatureMap::
initialize_cache_helix_stacking()
{
#ifdef USE_CACHE
  for (auto bp1: def_bps_)
  {
    auto ij1 = is_complementary_[bp1[0]][bp1[1]];
    for (auto bp2: def_bps_)
    {
      auto ij2 = is_complementary_[bp2[0]][bp2[1]];
      cache_helix_stacking_[ij1][ij2] = insert_key(s_helix_stacking+bp1+bp2);
    }
  }
#endif
}

size_t
FeatureMap::
find_helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
#ifdef USE_CACHE
  auto ij1 = is_complementary_[i1][j1];
  auto ij2 = is_complementary_[i2][j2];
  if (ij1>=0 && ij2>=0)
    return cache_helix_stacking_[ij1][ij2];
#endif
  return find_key(s_helix_stacking+i1+j1+i2+j2);
}

size_t
FeatureMap::
insert_helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2)
{
#ifdef USE_CACHE
  auto ij1 = is_complementary_[i1][j1];
  auto ij2 = is_complementary_[i2][j2];
  if (ij1>=0 && ij2>=0)
    return cache_helix_stacking_[ij1][ij2];
#endif
  return insert_key(s_helix_stacking+i1+j1+i2+j2);
}
#endif

#if PARAMS_HELIX_CLOSING
void
FeatureMap::
initialize_cache_helix_closing()
{
#ifdef USE_CACHE
  for (auto bp: def_bps_)
  {
    auto ij = is_complementary_[bp[0]][bp[1]];
    cache_helix_closing_[ij] = insert_key(s_helix_closing+bp);
  }
#endif
}

size_t
FeatureMap::
find_helix_closing(NUCL i, NUCL j) const
{
#ifdef USE_CACHE
  auto ij=is_complementary_[i][j];
  if (ij>=0)
    return cache_helix_closing_[ij];
#endif
  return find_key(s_helix_closing+i+j);
}

size_t
FeatureMap::
insert_helix_closing(NUCL i, NUCL j)
{
#ifdef USE_CACHE
  auto ij=is_complementary_[i][j];
  if (ij>=0)
    return cache_helix_closing_[ij];
#endif
  return insert_key(s_helix_closing+i+j);
}
#endif

#if PARAMS_MULTI_LENGTH
void
FeatureMap::
initialize_cache_multi_base()
{
#ifdef USE_CACHE
  cache_multi_base_ = insert_key(s_multi_base);
#endif
}

size_t
FeatureMap::
find_multi_base() const
{
#ifdef USE_CACHE
  return cache_multi_base_;
#else
  return find_key(s_multi_base);
#endif
}

size_t
FeatureMap::
insert_multi_base()
{
#ifdef USE_CACHE
  return cache_multi_base_;
#else
  return insert_key(s_multi_base);
#endif
}

void
FeatureMap::
initialize_cache_multi_unpaired()
{
#ifdef USE_CACHE
  cache_multi_unpaired_ = insert_key(s_multi_unpaired);
#endif
}

size_t
FeatureMap::
find_multi_unpaired() const
{
#ifdef USE_CACHE
  return cache_multi_unpaired_;
#else
  return find_key(s_multi_unpaired);
#endif
}

size_t
FeatureMap::
insert_multi_unpaired()
{
#ifdef USE_CACHE
  return cache_multi_unpaired_;
#else
  return insert_key(s_multi_unpaired);
#endif
}

void
FeatureMap::
initialize_cache_multi_paired()
{
#ifdef USE_CACHE
  cache_multi_paired_ = insert_key(s_multi_paired);
#endif
}

size_t
FeatureMap::
find_multi_paired() const
{
#ifdef USE_CACHE
  return cache_multi_paired_;
#else
  return find_key(s_multi_paired);
#endif
}

size_t
FeatureMap::
insert_multi_paired()
{
#ifdef USE_CACHE
  return cache_multi_paired_;
#else
  return insert_key(s_multi_paired);
#endif
}
#endif

#if PARAMS_DANGLE
void
FeatureMap::
initialize_cache_dangle_left()
{
#ifdef USE_CACHE
  for (auto bp: def_bps_)
  {
    auto ij = is_complementary_[bp[0]][bp[1]];
    for (auto i2: def_bases_)
    {
      auto ii2 = is_base_[i2];
      cache_dangle_left_[ij][ii2] = insert_key(s_dangle_left+bp+i2);
    }
  }
#endif
}

size_t
FeatureMap::
find_dangle_left(NUCL i1, NUCL j1, NUCL i2) const
{
#ifdef USE_CACHE
  auto ij1=is_complementary_[i1][j1];
  auto ii2=is_base_[i2];
  if (ij1>=0)
    return cache_dangle_left_[ij1][ii2];
#endif
  return find_key(s_dangle_left+i1+j1+i2);
}

size_t
FeatureMap::
insert_dangle_left(NUCL i1, NUCL j1, NUCL i2)
{
#ifdef USE_CACHE
  auto ij1=is_complementary_[i1][j1];
  auto ii2=is_base_[i2];
  if (ij1>=0)
    return cache_dangle_left_[ij1][ii2];
#endif
  return insert_key(s_dangle_left+i1+j1+i2);
}

void
FeatureMap::
initialize_cache_dangle_right()
{
#ifdef USE_CACHE
  for (auto bp: def_bps_)
  {
    auto ij = is_complementary_[bp[0]][bp[1]];
    for (auto j2: def_bases_)
    {
      auto jj2 = is_base_[j2];
      cache_dangle_right_[ij][jj2] = insert_key(s_dangle_right+bp+j2);
    }
  }
#endif
}

size_t
FeatureMap::
find_dangle_right(NUCL i1, NUCL j1, NUCL j2) const
{
#ifdef USE_CACHE
  auto ij1=is_complementary_[i1][j1];
  auto jj2=is_base_[j2];
  if (ij1>=0 && jj2>=0)
    return cache_dangle_right_[ij1][jj2];
#endif
    return find_key(s_dangle_right+i1+j1+j2);
}

size_t
FeatureMap::
insert_dangle_right(NUCL i1, NUCL j1, NUCL j2)
{
#ifdef USE_CACHE
  auto ij1=is_complementary_[i1][j1];
  auto jj2=is_base_[j2];
  if (ij1>=0 && jj2>=0)
    return cache_dangle_right_[ij1][jj2];
#endif
  return insert_key(s_dangle_right+i1+j1+j2);
}
#endif

#if PARAMS_EXTERNAL_LENGTH
void
FeatureMap::
initialize_cache_external_unpaired()
{
#ifdef USE_CACHE
  cache_external_unpaired_ = insert_key(s_external_unpaired);
#endif
}

size_t
FeatureMap::
find_external_unpaired() const
{
#ifdef USE_CACHE
  return cache_external_unpaired_;
#else
  return find_key(s_external_unpaired);
#endif
}

size_t
FeatureMap::
insert_external_unpaired()
{
#ifdef USE_CACHE
  return cache_external_unpaired_;
#else
  return insert_key(s_external_unpaired);
#endif
}

void
FeatureMap::
initialize_cache_external_paired()
{
#ifdef USE_CACHE
  cache_external_paired_ = insert_key(s_external_paired);
#endif
}

size_t
FeatureMap::
find_external_paired() const
{
#ifdef USE_CACHE
  return cache_external_paired_;
#else
  return find_key(s_external_paired);
#endif
}

size_t
FeatureMap::
insert_external_paired()
{
#ifdef USE_CACHE
  return cache_external_paired_;
#else
  return insert_key(s_external_paired);
#endif
}
#endif
