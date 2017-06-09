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

private:
  typedef cedar::da<int> trie_t;

public:
  ParameterHash();
  ~ParameterHash() { }
  ParameterHash(ParameterHash&& other);
  
  ParameterHash& operator=(ParameterHash&& other)
  {
    trie_ = std::move(other.trie_);
    values_ = std::move(other.values_);
    return *this;
  }

  void initialize();

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

  bool is_basepair_feature(const std::string& f) const;
  bool is_context_feature(const std::string& f) const;
  bool is_basepair_context_feature(const std::string& f) const;

  // read-only iterator
  class iterator : public std::iterator<std::forward_iterator_tag, ValueT >
  {
  public:
    iterator(const trie_t* trie, const std::vector<ValueT>* values, size_t from=0, size_t len=0)
      : trie_(const_cast<trie_t*>(trie)), values_(const_cast<std::vector<ValueT>*>(values)), from_(from), len_(len)
    { }

    iterator(trie_t* trie, std::vector<ValueT>* values, size_t from=0, size_t len=0)
      : trie_(trie), values_(values), from_(from), len_(len)
    { }

    iterator& begin()
    {
      i_ = trie_->begin(from_, len_);
      return *this;
    }

    iterator& end()
    {
      i_ = trie_t::CEDAR_NO_PATH;
      from_ = -1u;
      len_ = 0;
      return *this;
    }

    std::string key() const
    {
      std::string s(len_, ' ');
      trie_->suffix(&s[0], len_, from_);
      return s;
    }

    int index() const { return i_; }

    ValueT operator*() const
    {
      return (*values_)[i_];
    }

    ValueT& operator*()
    {
      return (*values_)[i_];
    }

    iterator& operator++()
    {
      i_ = trie_->next(from_, len_);
      if (i_==trie_t::CEDAR_NO_PATH) end();
      return *this;
    }

    iterator operator++(int)
    {
      iterator temp = *this;
      i_ = trie_->next(from_, len_);
      if (i_==trie_t::CEDAR_NO_PATH) end();
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
    std::vector<ValueT>* values_;
    size_t from_;
    size_t len_;
    int i_;
  };

  iterator begin()
  {
    return iterator(&trie_, &values_).begin();
  }

  iterator end()
  {
    return iterator(&trie_, &values_).end();
  }

  iterator cbegin() const
  {
    return iterator(&trie_, &values_).begin();
  }

  iterator cend() const
  {
    return iterator(&trie_, &values_).end();
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

  VI hairpin_nucleotides_cache(const std::vector<NUCL>& s, uint i, uint max_l) const;
  ValueT  hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l, const VI& pos) const;
  ValueT& hairpin_nucleotides(const std::vector<NUCL>& s, uint i, uint l, const VI& pos);

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
  ValueT  internal_nucleotides(const std::vector<NUCL>& s, const std::vector<NUCL>& t) const;
  ValueT& internal_nucleotides(const std::vector<NUCL>& s, const std::vector<NUCL>& t);

  VVI internal_nucleotides_cache(const std::vector<NUCL>& s, uint i, uint j, uint max_l, uint max_m) const;
  ValueT  internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m, const VVI& pos) const;
  ValueT& internal_nucleotides(const std::vector<NUCL>& s, uint i, uint l, uint j, uint m, const VVI& pos);

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
  static std::string def_bases;

  //std::unordered_map<std::string, ValueT> param_;
  trie_t trie_;
  std::vector<ValueT> values_;
  std::array<std::array<int, 256>, 256> is_complementary_;
  std::array<int, 256> is_base_;

  // cache
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
  int cache_internal_explicit_;
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
};

#endif  // PARAMETERHASH_HPP

// Local Variables:
// mode: C++
// End:

