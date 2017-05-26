#ifndef PARAMETERHASH_HPP
#define PARAMETERHASH_HPP

#include "Config.hpp"
#include "Utilities.hpp"
#include <string>
#include <unordered_map>
#include <iterator>
#include <array>
#include "cedar.h"

template < class ValueT >
class ParameterHash
{
public:
  //typedef typename ValueT ValueT;
  typedef cedar::da<ValueT> trie_t;

public:
  ParameterHash() { initialize(); }
  ~ParameterHash() { }

  ParameterHash(ParameterHash&& other)
    : param_(std::move(other.param_))
  {
    initialize();
  }

  ParameterHash& operator=(ParameterHash&& other)
  {
    param_ = std::move(other.param_);
    return *this;
  }

  void initialize();

  bool is_complementary(NUCL x, NUCL y) const;
  bool is_base(NUCL x) const;

  void LoadFromHash(const std::unordered_map<std::string, ValueT>& hash);
  void LoadDefaultComplementary();
  void LoadDefaultNonComplementary();
  void ReadFromFile(const std::string& filename);
  void WriteToFile(const std::string& filename, bool sort=true) const;

  ValueT  get_by_key(const std::string& key) const;
  ValueT& get_by_key(const std::string& key);

  bool is_basepair_feature(const std::string& f) const;
  bool is_context_feature(const std::string& f) const;
  bool is_basepair_context_feature(const std::string& f) const;

  // read-only iterator
  class iterator : public std::iterator<std::forward_iterator_tag, ValueT >
  {
  public:
    iterator(trie_t* trie, size_t from=0, size_t len=0)
      : trie_(trie), from_(from), len_(len)
    { }

    iterator(const trie_t* trie, size_t from=0, size_t len=0)
      : trie_(const_cast<trie_t*>(trie)), from_(from), len_(len)
    { }

    iterator& begin()
    {
      b_.i = trie_->begin(from_, len_);
      return *this;
    }

    iterator& end()
    {
      b_.i = trie_t::CEDAR_NO_PATH;
      from_ = -1u;
      len_ = 0;
      return *this;
    }

    std::string key() const
    {
      std::string s(len_+1, ' ');
      trie_->suffix(&s[0], len_, from_);
      return s;
    }

    ValueT operator*() const
    {
      return b_.x;
    }

    ValueT& operator*()
    {
      return trie_->update(key().c_str(), from_, len_, len_);
    }

    iterator& operator++()
    {
      b_.i = trie_->next(from_, len_);
      if (b_.i==trie_t::CEDAR_NO_PATH) end();
      return *this;
    }

    iterator operator++(int)
    {
      iterator temp = *this;
      b_.i = trie_->next(from_, len_);
      if (b_.i==trie_t::CEDAR_NO_PATH) end();
      return temp;
    }

    bool operator!=(const iterator& x) const
    {
      return from_!=x.from_ || len_!=x.len_;
    }

    bool operator==(const iterator& x) const
    {
      return from_==x.from_ && len_==x.len_;
    }

  private:
    trie_t* trie_;
    size_t from_;
    size_t len_;
    union { int i; typename trie_t::result_type x; } b_;
  };


  iterator begin()
  {
    return iterator(&param_).begin();
  }

  iterator end()
  {
    return iterator(&param_).end();
  }

  iterator cbegin() const
  {
    return iterator(&param_).begin();
  }

  iterator cend() const
  {
    return iterator(&param_).end();
  }

  //typename std::unordered_map<std::string, ValueT>::iterator begin() { return param_.begin(); }
  //typename std::unordered_map<std::string, ValueT>::const_iterator begin() const { return param_.begin(); }
  //typename std::unordered_map<std::string, ValueT>::iterator end() { return param_.end(); }
  //typename std::unordered_map<std::string, ValueT>::const_iterator end() const { return param_.end(); }
  //typename std::unordered_map<std::string, ValueT>::iterator erase(typename std::unordered_map<std::string, ValueT>::const_iterator pos) { return param_.erase(pos); }

  // access to parameters
#if PARAMS_BASE_PAIR
  ValueT  base_pair(NUCL i, NUCL j) const;
  ValueT& base_pair(NUCL i, NUCL j);
#endif
#if PARAMS_BASE_PAIR_DIST
  ValueT  base_pair_dist_at_least(uint l) const;
  ValueT& base_pair_dist_at_least(uint l);
#endif
#if PARAMS_TERMINAL_MISMATCH
  ValueT  terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const;
  ValueT& terminal_mismatch(NUCL i1, NUCL j1, NUCL i2, NUCL j2);  
#endif
#if PARAMS_HAIRPIN_LENGTH
  ValueT  hairpin_length_at_least(uint l) const;
  ValueT& hairpin_length_at_least(uint l);
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
  ValueT  hairpin_3_nucleotides(NUCL i1, NUCL i2, NUCL i3) const;
  ValueT& hairpin_3_nucleotides(NUCL i1, NUCL i2, NUCL i3);  
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
  ValueT  hairpin_4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const;
  ValueT& hairpin_4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4);  
#endif
#if PARAMS_HAIRPIN_5_NUCLEOTIDES
  ValueT  hairpin_5_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const;
  ValueT& hairpin_5_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5);  
#endif
#if PARAMS_HAIRPIN_6_NUCLEOTIDES
  ValueT  hairpin_6_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const;
  ValueT& hairpin_6_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6);  
#endif
#if PARAMS_HAIRPIN_7_NUCLEOTIDES
  ValueT  hairpin_7_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const;
  ValueT& hairpin_7_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7);  
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
  ValueT  bulge_0x1_nucleotides(NUCL i1) const;
  ValueT& bulge_0x1_nucleotides(NUCL i1);
  ValueT  bulge_1x0_nucleotides(NUCL i1) const;
  ValueT& bulge_1x0_nucleotides(NUCL i1);
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
  ValueT  bulge_0x2_nucleotides(NUCL i1, NUCL i2) const;
  ValueT& bulge_0x2_nucleotides(NUCL i1, NUCL i2);
  ValueT  bulge_2x0_nucleotides(NUCL i1, NUCL i2) const;
  ValueT& bulge_2x0_nucleotides(NUCL i1, NUCL i2);
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
  ValueT  bulge_0x3_nucleotides(NUCL i1, NUCL i2, NUCL i3) const;
  ValueT& bulge_0x3_nucleotides(NUCL i1, NUCL i2, NUCL i3);
  ValueT  bulge_3x0_nucleotides(NUCL i1, NUCL i2, NUCL i3) const;
  ValueT& bulge_3x0_nucleotides(NUCL i1, NUCL i2, NUCL i3);
#endif
#if PARAMS_BULGE_0x4_NUCLEOTIDES
  ValueT  bulge_0x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const;
  ValueT& bulge_0x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4);
  ValueT  bulge_4x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const;
  ValueT& bulge_4x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4);
#endif
#if PARAMS_BULGE_0x5_NUCLEOTIDES
  ValueT  bulge_0x5_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const;
  ValueT& bulge_0x5_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5);
  ValueT  bulge_5x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const;
  ValueT& bulge_5x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5);
#endif
#if PARAMS_BULGE_0x6_NUCLEOTIDES
  ValueT  bulge_0x6_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const;
  ValueT& bulge_0x6_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6);
  ValueT  bulge_6x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const;
  ValueT& bulge_6x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6);
#endif
#if PARAMS_BULGE_0x7_NUCLEOTIDES
  ValueT  bulge_0x7_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const;
  ValueT& bulge_0x7_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7);
  ValueT  bulge_7x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const;
  ValueT& bulge_7x0_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7);
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
  ValueT  internal_1x1_nucleotides(NUCL i1, NUCL i2) const;
  ValueT& internal_1x1_nucleotides(NUCL i1, NUCL i2);
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
  ValueT  internal_1x2_nucleotides(NUCL i1, NUCL i2, NUCL i3) const;
  ValueT& internal_1x2_nucleotides(NUCL i1, NUCL i2, NUCL i3);
  ValueT  internal_2x1_nucleotides(NUCL i1, NUCL i2, NUCL i3) const;
  ValueT& internal_2x1_nucleotides(NUCL i1, NUCL i2, NUCL i3);
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
  ValueT  internal_2x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const;
  ValueT& internal_2x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4);
#endif
#if PARAMS_INTERNAL_1x3_NUCLEOTIDES
  ValueT  internal_1x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const;
  ValueT& internal_1x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4);
  ValueT  internal_3x1_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4) const;
  ValueT& internal_3x1_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4);
#endif
#if PARAMS_INTERNAL_2x3_NUCLEOTIDES
  ValueT  internal_2x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const;
  ValueT& internal_2x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5);
  ValueT  internal_3x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const;
  ValueT& internal_3x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5);
#endif
#if PARAMS_INTERNAL_3x3_NUCLEOTIDES
  ValueT  internal_3x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const;
  ValueT& internal_3x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6);
#endif

#if PARAMS_INTERNAL_1x4_NUCLEOTIDES
  ValueT  internal_1x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const;
  ValueT& internal_1x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5);
  ValueT  internal_4x1_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5) const;
  ValueT& internal_4x1_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5);
#endif
#if PARAMS_INTERNAL_2x4_NUCLEOTIDES
  ValueT  internal_2x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const;
  ValueT& internal_2x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6);
  ValueT  internal_4x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6) const;
  ValueT& internal_4x2_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6);
#endif
#if PARAMS_INTERNAL_3x4_NUCLEOTIDES
  ValueT  internal_3x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const;
  ValueT& internal_3x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7);
  ValueT  internal_4x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7) const;
  ValueT& internal_4x3_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7);
#endif
#if PARAMS_INTERNAL_4x4_NUCLEOTIDES
  ValueT  internal_4x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7, NUCL i8) const;
  ValueT& internal_4x4_nucleotides(NUCL i1, NUCL i2, NUCL i3, NUCL i4, NUCL i5, NUCL i6, NUCL i7, NUCL i8);
#endif
#if PARAMS_HELIX_STACKING
  ValueT  helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2) const;
  ValueT& helix_stacking(NUCL i1, NUCL j1, NUCL i2, NUCL j2);
#endif
#if PARAMS_HELIX_CLOSING
  ValueT  helix_closing(NUCL i, NUCL j) const;
  ValueT& helix_closing(NUCL i, NUCL j);
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
  ValueT  dangle_left(NUCL i1, NUCL j1, NUCL i2) const;
  ValueT& dangle_left(NUCL i1, NUCL j1, NUCL i2);
  ValueT  dangle_right(NUCL i1, NUCL j1, NUCL j2) const;
  ValueT& dangle_right(NUCL i1, NUCL j1, NUCL j2);
#endif
#if PARAMS_EXTERNAL_LENGTH
  ValueT  external_unpaired() const;
  ValueT& external_unpaired();
  ValueT  external_paired() const;
  ValueT& external_paired();
#endif

private:
  //std::unordered_map<std::string, ValueT> param_;
  trie_t param_;
  std::array<std::array<int, 256>, 256> is_complementary_;
  std::array<int, 256> is_base_;
};

#endif  // PARAMETERHASH_HPP

// Local Variables:
// mode: C++
// End:

