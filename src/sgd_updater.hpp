#ifndef __INC_SGD_UPDATER_HPP__
#define __INC_SGD_UPDATER_HPP__

#include <string>
#include <unordered_map>

class SGDUpdater
{
public:
  SGDUpdater(double eta);
  
  double update(const std::string& fname, double& w, double grad);

private:
  double eta_;
};

class AdaGradRDAUpdater
{
public:
  AdaGradRDAUpdater(double eta,  double lambda, double eps=1e-8);

  double update(const std::string& fname, double& w, double grad);
  void proceed_time() { ++t_; }

private:
  double eta_;
  double lambda_;
  double eps_;
  unsigned int t_;
  std::unordered_map<std::string,double> sum_grad_;
  std::unordered_map<std::string,double> sum_squared_grad_;
};

#endif //  __INC_SGD_UPDATER_HPP__
