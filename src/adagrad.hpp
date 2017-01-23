#ifndef __INC_SGD_UPDATER_HPP__
#define __INC_SGD_UPDATER_HPP__

#include <string>
#include <unordered_map>
#include "ParameterHash.hpp"

class AdaGradRDAUpdater
{
public:
  AdaGradRDAUpdater(double eta,  double lambda, double eps=1e-8);

  void update(const std::string& fname, double& w, double grad);
  void regularize(const std::string& fname, double& w) const;
  void proceed_timestamp() { ++t_; }

  void read_from_file(const std::string& filename);
  void write_to_file(const std::string& filename, const ParameterHash<double>* pm, bool sort=true) const;

  std::unordered_map<std::string,double>::const_iterator begin() const { return sum_grad_.begin(); }
  std::unordered_map<std::string,double>::const_iterator end() const { return sum_grad_.end(); }

private:
  double eta_;
  double lambda_;
  double eps_;
  unsigned int t_;
  std::unordered_map<std::string,double> sum_grad_;
  std::unordered_map<std::string,double> sum_squared_grad_;
};

class AdaGradFobosUpdater
{
public:
  AdaGradFobosUpdater(double eta,  double lambda, double eps=1e-8);

  void update(const std::string& fname, double& w, double grad, double weight);
  void regularize(const std::string& fname, double& w, double weight) const;
  void proceed_timestamp() { }

  void read_from_file(const std::string& filename);
  void write_to_file(const std::string& filename, const ParameterHash<double>* pm, bool sort=true) const;

private:
  double eta_;
  double lambda_;
  double eps_;
  std::unordered_map<std::string,double> sum_squared_grad_;
};

#endif //  __INC_SGD_UPDATER_HPP__
