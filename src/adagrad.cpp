#include "../config.h"
#include <cmath>
#include "adagrad.hpp"

inline
double
clip(double w, double c)
{
  if (w>=0.0)
    return w>c ? w-c : 0.0;
  else
    return -clip(-w, c);
}

AdaGradRDAUpdater::
AdaGradRDAUpdater(double eta, double lambda, double eps)
  : eta_(eta), lambda_(lambda), eps_(eps), t_(1), sum_grad_(), sum_squared_grad_()
{
}

void
AdaGradRDAUpdater::
update(const std::string& fname, double& w, double grad)
{
  auto u = sum_grad_.insert(std::make_pair(fname, 0.0));
  u.first->second += grad;
  auto g = sum_squared_grad_.insert(std::make_pair(fname, eps_));
  g.first->second += grad*grad;
  if (!u.second) regularize(fname, w);
}

void
AdaGradRDAUpdater::
regularize(const std::string& fname, double& w) const
{
  auto u = sum_grad_.find(fname);
  auto g = sum_squared_grad_.find(fname);
  w = - eta_ * t_ / std::sqrt(g->second) * clip(u->second/t_, lambda_);
}


AdaGradFobosUpdater::
AdaGradFobosUpdater(double eta, double lambda, double eps)
  : eta_(eta), lambda_(lambda), eps_(eps), sum_squared_grad_()
{
}

void
AdaGradFobosUpdater::
update(const std::string& fname, double& w, double grad)
{
  auto g = sum_squared_grad_.insert(std::make_pair(fname, eps_));
  g.first->second += grad*grad;
  w -= eta_ / std::sqrt(g.first->second) * grad;
}

void
AdaGradFobosUpdater::
regularize(const std::string& fname, double& w) const
{
  w = clip(w, lambda_);
}
