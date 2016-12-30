#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <cerrno>
#include <cstdio>
#include <cctype>
#include "ParameterHash.hpp"
#include "Utilities.hpp"

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
LoadFromHash(std::unordered_map<std::string, ValueT> hash)
{
  param_ = std::move(hash);
}

template < class ValueT >
void
ParameterHash<ValueT>::
ReadFromFile(const std::string& filename)
{
  param_.clear();
  std::ifstream is(filename.c_str());
  if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);
  std::string k;
  ValueT v;
  while (is >> k >> v)
  {
    if (v!=0.0)
      param_.insert(std::make_pair(k, v));
  }
}

template < class ValueT >
void
ParameterHash<ValueT>::
WriteToFile(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
  if (!os) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);
  
  for (const auto& e : param_)
  {
    if (e.second!=0.0)
      os << e.first << " " << e.second << std::endl;
  }
}

template < class ValueT >
inline
ValueT
ParameterHash<ValueT>::
get_by_key(const std::string& key) const
{
  auto itr = param_.find(key);
  return itr==param_.end() ? static_cast<ValueT>(0) : itr->second;
}

template < class ValueT >
inline
ValueT&
ParameterHash<ValueT>::
get_by_key(const std::string& key)
{
  auto ret = param_.insert(std::make_pair(key, static_cast<ValueT>(0)));
  return ret.first->second;
}

#if PARAMS_BASE_PAIR
static const char* format_base_pair = "base_pair_%c%c";

template < class ValueT >
inline
ValueT
ParameterHash<ValueT>::
base_pair(NUCL i, NUCL j) const
{
  auto v = std::minmax(i, j);
  return get_by_key(string_format(format_base_pair, v.first, v.second));
}

template < class ValueT >
inline
ValueT&
ParameterHash<ValueT>::
base_pair(NUCL i, NUCL j)
{
  auto v = std::minmax(i, j);
  return get_by_key(string_format(format_base_pair, v.first, v.second));
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

#if PARAMS_HAIRPIN_3_NUCLEOTIDES
static const char* format_hairpin_3_nucleotides = "hairpin_3_nucleotides_%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_3_nucleotides(NUCL i1, NUCL i2, NUCL i3) const
{
  return get_by_key(string_format(format_hairpin_3_nucleotides, i1, i2, i3));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_3_nucleotides(NUCL i1, NUCL i2, NUCL i3)
{
  return get_by_key(string_format(format_hairpin_3_nucleotides, i1, i2, i3));
}
#endif

#if PARAMS_HAIRPIN_4_NUCLEOTIDES
static const char* format_hairpin_4_nucleotides = "hairpin_4_nucleotides_%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const
{
  return get_by_key(string_format(format_hairpin_4_nucleotides, i1, i2, i3, i4));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4)
{
  return get_by_key(string_format(format_hairpin_4_nucleotides, i1, i2, i3, i4));
}
#endif

#if PARAMS_HAIRPIN_5_NUCLEOTIDES
static const char* format_hairpin_5_nucleotides = "hairpin_5_nucleotides_%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_5_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const
{
  return get_by_key(string_format(format_hairpin_5_nucleotides, i1, i2, i3, i4, i5));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_5_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5)
{
  return get_by_key(string_format(format_hairpin_5_nucleotides, i1, i2, i3, i4, i5));
}
#endif

#if PARAMS_HAIRPIN_6_NUCLEOTIDES
static const char* format_hairpin_6_nucleotides = "hairpin_6_nucleotides_%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_6_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const
{
  return get_by_key(string_format(format_hairpin_6_nucleotides, i1, i2, i3, i4, i5, i6));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_6_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6)
{
  return get_by_key(string_format(format_hairpin_6_nucleotides, i1, i2, i3, i4, i5, i6));
}
#endif

#if PARAMS_HAIRPIN_7_NUCLEOTIDES
static const char* format_hairpin_7_nucleotides = "hairpin_7_nucleotides_%c%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
hairpin_7_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const
{
  return get_by_key(string_format(format_hairpin_7_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
hairpin_7_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7)
{
  return get_by_key(string_format(format_hairpin_7_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}
#endif

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
  auto v = std::minmax(i, j);
  return get_by_key(string_format(format_internal_explicit, v.first, v.second));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_explicit(uint i, uint j)
{
  auto v = std::minmax(i, j);
  return get_by_key(string_format(format_internal_explicit, v.first, v.second));
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

#if PARAMS_BULGE_0x1_NUCLEOTIDES
static const char* format_bulge_0x1_nucleotides = "bulge_0x1_nucleotides_%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_0x1_nucleotides(NUCL i1) const
{
  return get_by_key(string_format(format_bulge_0x1_nucleotides, i1));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_0x1_nucleotides(NUCL i1)
{
  return get_by_key(string_format(format_bulge_0x1_nucleotides, i1));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_1x0_nucleotides(NUCL i1) const
{
  return get_by_key(string_format(format_bulge_0x1_nucleotides, i1));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_1x0_nucleotides(NUCL i1)
{
  return get_by_key(string_format(format_bulge_0x1_nucleotides, i1));
}
#endif

#if PARAMS_BULGE_0x2_NUCLEOTIDES
static const char* format_bulge_0x2_nucleotides = "bulge_0x2_nucleotides_%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_0x2_nucleotides(NUCL i1, NUCL i2) const
{
  return get_by_key(string_format(format_bulge_0x2_nucleotides, i1, i2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_0x2_nucleotides(NUCL i1, NUCL i2)
{
  return get_by_key(string_format(format_bulge_0x2_nucleotides, i1, i2));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_2x0_nucleotides(NUCL i1, NUCL i2) const
{
  return get_by_key(string_format(format_bulge_0x2_nucleotides, i1, i2));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_2x0_nucleotides(NUCL i1, NUCL i2)
{
  return get_by_key(string_format(format_bulge_0x2_nucleotides, i1, i2));
}
#endif

#if PARAMS_BULGE_0x3_NUCLEOTIDES
static const char* format_bulge_0x3_nucleotides = "bulge_0x3_nucleotides_%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_0x3_nucleotides(NUCL i1, NUCL i2, NUCL i3) const
{
  return get_by_key(string_format(format_bulge_0x3_nucleotides, i1, i2, i3));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_0x3_nucleotides(NUCL i1, NUCL i2, NUCL i3)
{
  return get_by_key(string_format(format_bulge_0x3_nucleotides, i1, i2, i3));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_3x0_nucleotides(NUCL i1, NUCL i2, NUCL i3) const
{
  return get_by_key(string_format(format_bulge_0x3_nucleotides, i1, i2, i3));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_3x0_nucleotides(NUCL i1, NUCL i2, NUCL i3)
{
  return get_by_key(string_format(format_bulge_0x3_nucleotides, i1, i2, i3));
}
#endif

#if PARAMS_BULGE_0x4_NUCLEOTIDES
static const char* format_bulge_0x4_nucleotides = "bulge_0x4_nucleotides_%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_0x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const
{
  return get_by_key(string_format(format_bulge_0x4_nucleotides, i1, i2, i3, i4));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_0x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4)
{
  return get_by_key(string_format(format_bulge_0x4_nucleotides, i1, i2, i3, i4));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_4x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const
{
  return get_by_key(string_format(format_bulge_0x4_nucleotides, i1, i2, i3, i4));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_4x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4)
{
  return get_by_key(string_format(format_bulge_0x4_nucleotides, i1, i2, i3, i4));
}
#endif

#if PARAMS_BULGE_0x5_NUCLEOTIDES
static const char* format_bulge_0x5_nucleotides = "bulge_0x5_nucleotides_%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_0x5_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const
{
  return get_by_key(string_format(format_bulge_0x5_nucleotides, i1, i2, i3, i4, i5));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_0x5_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5)
{
  return get_by_key(string_format(format_bulge_0x5_nucleotides, i1, i2, i3, i4, i5));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_5x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const
{
  return get_by_key(string_format(format_bulge_0x5_nucleotides, i1, i2, i3, i4, i5));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_5x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5)
{
  return get_by_key(string_format(format_bulge_0x5_nucleotides, i1, i2, i3, i4, i5));
}
#endif

#if PARAMS_BULGE_0x6_NUCLEOTIDES
static const char* format_bulge_0x6_nucleotides = "bulge_0x6_nucleotides_%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_0x6_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const
{
  return get_by_key(string_format(format_bulge_0x6_nucleotides, i1, i2, i3, i4, i5, i6));
}
  

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_0x6_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6)
{
  return get_by_key(string_format(format_bulge_0x6_nucleotides, i1, i2, i3, i4, i5, i6));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_6x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const
{
  return get_by_key(string_format(format_bulge_0x6_nucleotides, i1, i2, i3, i4, i5, i6));
}
 
template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_6x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6)
{
  return get_by_key(string_format(format_bulge_0x6_nucleotides, i1, i2, i3, i4, i5, i6));
}
#endif

#if PARAMS_BULGE_0x7_NUCLEOTIDES
static const char* format_bulge_0x7_nucleotides = "bulge_0x7_nucleotides_%c%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_0x7_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const
{
  return get_by_key(string_format(format_bulge_0x7_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_0x7_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7)
{
  return get_by_key(string_format(format_bulge_0x7_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
bulge_7x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const
{
  return get_by_key(string_format(format_bulge_0x7_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
bulge_7x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7)
{
  return get_by_key(string_format(format_bulge_0x7_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}
#endif

#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
static const char* format_internal_1x1_nucleotides = "internal_1x1_nucleotides_%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_1x1_nucleotides(NUCL i1, NUCL i2) const
{
  auto s1 = string_format(format_internal_1x1_nucleotides, i1, i2);
  auto s2 = string_format(format_internal_1x1_nucleotides, i2, i1);
  return get_by_key(s1<s2 ? s1 : s2);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_1x1_nucleotides(NUCL i1, NUCL i2)
{
  auto s1 = string_format(format_internal_1x1_nucleotides, i1, i2);
  auto s2 = string_format(format_internal_1x1_nucleotides, i2, i1);
  return get_by_key(s1<s2 ? s1 : s2);
}
#endif

#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
static const char* format_internal_1x2_nucleotides = "internal_1x2_nucleotides_%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_1x2_nucleotides(NUCL i1, NUCL i2, NUCL i3) const
{
  return get_by_key(string_format(format_internal_1x2_nucleotides, i1, i2, i3));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_1x2_nucleotides(NUCL i1, NUCL i2, NUCL i3)
{
  return get_by_key(string_format(format_internal_1x2_nucleotides, i1, i2, i3));
}

static const char* format_internal_2x1_nucleotides = "internal_2x1_nucleotides_%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_2x1_nucleotides(NUCL i1, NUCL i2, NUCL i3) const
{
  return get_by_key(string_format(format_internal_2x1_nucleotides, i1, i2, i3));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_2x1_nucleotides(NUCL i1, NUCL i2, NUCL i3)
{
  return get_by_key(string_format(format_internal_2x1_nucleotides, i1, i2, i3));
}
#endif

#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
static const char* format_internal_2x2_nucleotides = "internal_2x2_nucleotides_%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_2x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const
{
  auto s1 = string_format(format_internal_2x2_nucleotides, i1, i2, i3, i4);
  auto s2 = string_format(format_internal_2x2_nucleotides, i3, i4, i1, i2);
  return get_by_key(s1<s2 ? s1 : s2);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_2x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4)
{
  auto s1 = string_format(format_internal_2x2_nucleotides, i1, i2, i3, i4);
  auto s2 = string_format(format_internal_2x2_nucleotides, i3, i4, i1, i2);
  return get_by_key(s1<s2 ? s1 : s2);
}
#endif

#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
static const char* format_internal_1x3_nucleotides = "internal_1x3_nucleotides_%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_1x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const
{
  return get_by_key(string_format(format_internal_1x3_nucleotides, i1, i2, i3, i4));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_1x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4)
{
  return get_by_key(string_format(format_internal_1x3_nucleotides, i1, i2, i3, i4));
}

static const char* format_internal_3x1_nucleotides = "internal_3x1_nucleotides_%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_3x1_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const
{
  return get_by_key(string_format(format_internal_3x1_nucleotides, i1, i2, i3, i4));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_3x1_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4)
{
  return get_by_key(string_format(format_internal_3x1_nucleotides, i1, i2, i3, i4));
}
#endif

#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
static const char* format_internal_2x3_nucleotides = "internal_2x3_nucleotides_%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_2x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const
{
  return get_by_key(string_format(format_internal_2x3_nucleotides, i1, i2, i3, i4, i5));
}  

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_2x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5)
{
  return get_by_key(string_format(format_internal_2x3_nucleotides, i1, i2, i3, i4, i5));
}  

static const char* format_internal_3x2_nucleotides = "internal_3x2_nucleotides_%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_3x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const
{
  return get_by_key(string_format(format_internal_3x2_nucleotides, i1, i2, i3, i4, i5));
}  

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_3x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5)
{
  return get_by_key(string_format(format_internal_3x2_nucleotides, i1, i2, i3, i4, i5));
}  
#endif

#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
static const char* format_internal_3x3_nucleotides = "internal_3x3_nucleotides_%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_3x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const
{
  auto s1 = string_format(format_internal_3x3_nucleotides, i1, i2, i3, i4, i5, i6);
  auto s2 = string_format(format_internal_3x3_nucleotides, i4, i5, i6, i1, i2, i3);
  return get_by_key(s1<s2 ? s1 : s2);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_3x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6)
{
  auto s1 = string_format(format_internal_3x3_nucleotides, i1, i2, i3, i4, i5, i6);
  auto s2 = string_format(format_internal_3x3_nucleotides, i4, i5, i6, i1, i2, i3);
  return get_by_key(s1<s2 ? s1 : s2);
}
#endif

#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
static const char* format_internal_1x4_nucleotides = "internal_1x4_nucleotides_%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_1x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const
{
  return get_by_key(string_format(format_internal_1x4_nucleotides, i1, i2, i3, i4, i5));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_1x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5)
{
  return get_by_key(string_format(format_internal_1x4_nucleotides, i1, i2, i3, i4, i5));
}

static const char* format_internal_4x1_nucleotides = "internal_4x1_nucleotides_%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_4x1_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const
{
  return get_by_key(string_format(format_internal_4x1_nucleotides, i1, i2, i3, i4, i5));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_4x1_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5)
{
  return get_by_key(string_format(format_internal_4x1_nucleotides, i1, i2, i3, i4, i5));
}
#endif

#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
static const char* format_internal_2x4_nucleotides = "internal_2x4_nucleotides_%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_2x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const
{
  return get_by_key(string_format(format_internal_2x4_nucleotides, i1, i2, i3, i4, i5, i6));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_2x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6)
{
  return get_by_key(string_format(format_internal_2x4_nucleotides, i1, i2, i3, i4, i5, i6));
}

static const char* format_internal_4x2_nucleotides = "internal_4x2_nucleotides_%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_4x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const
{
  return get_by_key(string_format(format_internal_4x2_nucleotides, i1, i2, i3, i4, i5, i6));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_4x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6)
{
  return get_by_key(string_format(format_internal_4x2_nucleotides, i1, i2, i3, i4, i5, i6));
}
#endif

#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
static const char* format_internal_3x4_nucleotides = "internal_3x4_nucleotides_%c%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_3x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const
{
  return get_by_key(string_format(format_internal_3x4_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_3x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7)
{
  return get_by_key(string_format(format_internal_3x4_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}

static const char* format_internal_4x3_nucleotides = "internal_4x3_nucleotides_%c%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_4x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const
{
  return get_by_key(string_format(format_internal_4x3_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_4x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7)
{
  return get_by_key(string_format(format_internal_4x3_nucleotides, i1, i2, i3, i4, i5, i6, i7));
}
#endif

#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
static const char* format_internal_4x4_nucleotides = "internal_4x4_nucleotides_%c%c%c%c%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
internal_4x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7, NUCL i8) const
{
  auto s1 = string_format(format_internal_4x4_nucleotides, i1, i2, i3, i4, i5, i6, i7, i8);
  auto s2 = string_format(format_internal_4x4_nucleotides, i5, i6, i7, i8, i1, i2, i3, i4);
  return get_by_key(s1<s2 ? s1 : s2);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
internal_4x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7, NUCL i8)
{
  auto s1 = string_format(format_internal_4x4_nucleotides, i1, i2, i3, i4, i5, i6, i7, i8);
  auto s2 = string_format(format_internal_4x4_nucleotides, i5, i6, i7, i8, i1, i2, i3, i4);
  return get_by_key(s1<s2 ? s1 : s2);
}
#endif

#if PARAMS_HELIX_STACKING
static const char* format_helix_stacking = "helix_stacking_%c%c%c%c";

template <class ValueT>
inline
ValueT
ParameterHash<ValueT>::
helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const
{
  auto s1 = string_format(format_helix_stacking, i1, j1, i2, j2);
  auto s2 = string_format(format_helix_stacking, j2, i2, j1, i1);
  return get_by_key(s1<s2 ? s1 : s2);
}

template <class ValueT>
inline
ValueT&
ParameterHash<ValueT>::
helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2)
{
  auto s1 = string_format(format_helix_stacking, i1, j1, i2, j2);
  auto s2 = string_format(format_helix_stacking, j2, i2, j1, i1);
  return get_by_key(s1<s2 ? s1 : s2);
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

template
class ParameterHash<double>;

template
class ParameterHash<uint>;
