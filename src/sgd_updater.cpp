#include "../config.h"
#include <cmath>
#include "sgd_updater.hpp"

SGDUpdater::
SGDUpdater(double eta)
  : eta_(eta)
{
}

double
SGDUpdater::
update(const std::string& fname, double& w, double grad)
{
  auto v = eta_ * grad;
  w -= v;
  return v;
}

AdaGradRDAUpdater::
AdaGradRDAUpdater(double eta, double lambda, double eps)
  : eta_(eta), lambda_(lambda), eps_(eps), t_(1), sum_grad_(), sum_squared_grad_()
{
}

double
AdaGradRDAUpdater::
update(const std::string& fname, double& w, double grad)
{
  auto u = sum_grad_.insert(std::make_pair(fname, 0.0));
  u.first->second += grad;
  auto g = sum_squared_grad_.insert(std::make_pair(fname, eps_));
  g.first->second += grad*grad;
  if (u.first->second>0.0)
  {
    auto v = u.first->second/t_ - lambda_;
    w = v<=0.0 ? 0.0 : - eta_ * t_ / std::sqrt(g.first->second) * v;
    return v-w; 
  }
  else if (u.first->second<0.0)
  {
    auto v = -u.first->second/t_ - lambda_;
    w = v<=0.0 ? 0.0 : eta_ * t_ / std::sqrt(g.first->second) * v;
    return v-w; 
  }
  return 0.0;
}
