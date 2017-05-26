#ifndef __INC_SGD_UPDATER_HPP__
#define __INC_SGD_UPDATER_HPP__

#include <string>
#include <unordered_map>
#include "Config.hpp"
#include "ParameterHash.hpp"

class AdaGradRDAUpdater
{
public:
  AdaGradRDAUpdater(float eta,  float lambda, float eps=1e-8);

  void update(const std::string& fname, param_value_type& w, param_value_type grad, param_value_type weight);
  void regularize(const std::string& fname, param_value_type& w, param_value_type weight) const;
  void proceed_timestamp() { ++t_; }

  void read_from_file(const std::string& filename);
  void write_to_file(const std::string& filename, const ParameterHash<param_value_type>* pm, bool sort=true) const;

  std::unordered_map<std::string,param_value_type>::const_iterator begin() const { return sum_grad_.begin(); }
  std::unordered_map<std::string,param_value_type>::const_iterator end() const { return sum_grad_.end(); }

private:
  float eta_;
  float lambda_;
  float eps_;
  unsigned int t_;
  std::unordered_map<std::string,param_value_type> sum_grad_;
  std::unordered_map<std::string,param_value_type> sum_squared_grad_;
};

class AdaGradFobosUpdater
{
public:
  AdaGradFobosUpdater(float eta,  float lambda, float eps=1e-8);


  void update(const std::string& fname, param_value_type& w, param_value_type grad, param_value_type weight);
  void regularize(const std::string& fname, param_value_type& w, param_value_type weight) const;
  void proceed_timestamp() { }

  void read_from_file(const std::string& filename);
  void write_to_file(const std::string& filename, const ParameterHash<param_value_type>* pm, bool sort=true) const;

private:
  float eta_;
  float lambda_;
  float eps_;
  std::unordered_map<std::string,param_value_type> sum_squared_grad_;
};

#endif //  __INC_SGD_UPDATER_HPP__
