#include "../config.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cerrno>
#include <cstdio>
#include <cctype>
#include <cassert>
#include "adagrad.hpp"

template < class T >
inline
T
clip(T w, T c)
{
  if (w>=0.0)
    return w>c ? w-c : 0.0;
  else
    return -clip(-w, c);
}

AdaGradRDAUpdater::
AdaGradRDAUpdater(float eta, float lambda, float eps)
  : eta_(eta), lambda_(lambda), eps_(eps), t_(1), sum_grad_(), sum_squared_grad_()
{
}

void
AdaGradRDAUpdater::
update(const std::string& fname, param_value_type& w, param_value_type grad, param_value_type weight)
{
  auto u = sum_grad_.insert(std::make_pair(fname, 0.0));
  u.first->second += grad;
  auto g = sum_squared_grad_.insert(std::make_pair(fname, eps_));
  g.first->second += grad*grad;
  if (!u.second) regularize(fname, w, weight);
}

void
AdaGradRDAUpdater::
regularize(const std::string& fname, param_value_type& w, param_value_type weight) const
{
  auto u = sum_grad_.find(fname);
  auto g = sum_squared_grad_.find(fname);
  w = - weight * eta_ * t_ / std::sqrt(g->second) * clip(u->second/t_, static_cast<param_value_type>(lambda_));
}

void
AdaGradRDAUpdater::
read_from_file(const std::string& filename)
{
  sum_grad_.clear();
  sum_squared_grad_.clear();
  std::ifstream is(filename.c_str());
  if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  is >> t_ >> eta_ >> lambda_ >> eps_;

  std::string fname;
  double p, s1, s2;
  while (is >> fname >> p >> s1 >> s2)
  {
    sum_grad_[fname] = s1;
    sum_squared_grad_[fname] = s2;
  }
}

void
AdaGradRDAUpdater::
write_to_file(const std::string& filename, const ParameterHash<param_value_type>* pm, bool sort) const
{
  std::ofstream os(filename.c_str());
  if (!os) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  os << t_ << " " << eta_ << " " << lambda_ << " " << eps_ << std::endl;

  if (sort)
  {
    std::vector<std::string> keys;
    for (const auto& e : sum_grad_)
      keys.emplace_back(e.first);
    std::sort(keys.begin(), keys.end());
    for (const auto& k : keys)
    {
      os << k << " " <<  pm->get_by_key(k) << " " 
         << sum_grad_.find(k)->second << " " 
         << sum_squared_grad_.find(k)->second << std::endl;
    }
  }
  else
  {
    for (const auto& e : sum_grad_)
    {
      os << e.first << " " << pm->get_by_key(e.first)  << " " << e.second << " " 
         << sum_squared_grad_.find(e.first)->second << std::endl;
    }
  }
}



AdaGradFobosUpdater::
AdaGradFobosUpdater(float eta, float lambda, float eps)
  : eta_(eta), lambda_(lambda), eps_(eps), sum_squared_grad_()
{
}

void
AdaGradFobosUpdater::
update(const std::string& fname, param_value_type& w, param_value_type grad, param_value_type weight)
{
  auto g = sum_squared_grad_.insert(std::make_pair(fname, eps_));
  g.first->second += grad*grad;
  std::cout << "  " << fname << ": w=" << w << ", g=" << grad << ", g2s=" << g.first->second;
  w -= weight * eta_ / std::sqrt(g.first->second) * grad;
  std::cout << ", update=" << weight * eta_ / std::sqrt(g.first->second) * grad << ", w_new=" << w << std::endl;
}

void
AdaGradFobosUpdater::
regularize(const std::string& fname, param_value_type& w, param_value_type weight) const
{
  auto g = sum_squared_grad_.find(fname);
  if (g!=sum_squared_grad_.end() && g->second>0.0) 
    w = clip(w, weight * eta_ / std::sqrt(g->second) * lambda_);
}

void
AdaGradFobosUpdater::
read_from_file(const std::string& filename)
{
  sum_squared_grad_.clear();
  std::ifstream is(filename.c_str());
  if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  if (std::isdigit(is.peek())) 
  {
    uint t;
    is >> t >> eta_ >> lambda_ >> eps_;

    std::string fname;
    double p, s1, s2;
    while (is >> fname >> p >> s1 >> s2)
      if (p!=0.0 || s1!=0.0 || s2!=0.0)
        sum_squared_grad_[fname] = s2;
  }
  else
  {
    std::string fname;
    double p;
    while (is >> fname >> p)
      if (p!=0.0)
        sum_squared_grad_[fname] = 0.0;
  }
}

void
AdaGradFobosUpdater::
write_to_file(const std::string& filename, const ParameterHash<param_value_type>* pm, bool sort) const
{
  std::ofstream os(filename.c_str());
  if (!os) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  os << 0 << " " << eta_ << " " << lambda_ << " " << eps_ << std::endl;

  if (sort)
  {
    std::vector<std::string> keys;
    for (const auto& e : sum_squared_grad_)
      keys.emplace_back(e.first);
    std::sort(keys.begin(), keys.end());
    for (const auto& k : keys)
    {
      os << k << " " <<  pm->get_by_key(k) << " " 
         << 0.0 << " " 
         << sum_squared_grad_.find(k)->second << std::endl;
    }
  }
  else
  {
    for (const auto& e : sum_squared_grad_)
    {
      os << e.first << " " << pm->get_by_key(e.first) << " " << 0.0 << " " 
         << e.second << std::endl;
    }
  }
}
