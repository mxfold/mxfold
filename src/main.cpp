#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <string>
#include <stdexcept>
#include "../config.h"
#include "cmdline.h"
#include "Config.hpp"
#include "Utilities.hpp"
#include "InferenceEngine.hpp"
#include "ParameterHash.hpp"
#include "SStruct.hpp"
#include "sgd_updater.hpp"

extern std::unordered_map<std::string, double> default_params_complementary;
extern std::unordered_map<std::string, double> default_params_noncomplementary;

class NGSfold
{
public:
  NGSfold() : train_mode_(false), mea_(false), gce_(false) { }

  NGSfold& parse_options(int& argc, char**& argv);

  int run()
  {
    if (train_mode_)
      return train();
    else
      return predict();
  }

private:
  int train();
  int predict();

private:
  bool train_mode_;
  bool noncomplementary_;
  bool output_bpseq_;
  std::string out_file_;
  std::string param_file_;
  std::vector<std::string> data_str_list_;
  std::vector<std::string> data_unpaired_list_;
  std::vector<std::string> data_paired_list_;  
  bool mea_;
  bool gce_;
  std::vector<float> gamma_;
  float pos_w_;
  float neg_w_;
  float lambda_;
  float eta0_;
  std::vector<std::string> args_;
};

NGSfold& 
NGSfold::parse_options(int& argc, char**& argv)
{
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info)!=0) exit(1);

  if (args_info.param_given)
    param_file_ = args_info.param_arg;

  if (args_info.train_given)
  {
    train_mode_ = true;
    out_file_ = args_info.train_arg;
  }

  if (args_info.mea_given)
  {
    mea_ = true;
    gamma_.resize(args_info.mea_given);
    for (uint i=0; i!=args_info.mea_given; ++i)
      gamma_[i] = args_info.mea_arg[i];
  }

  if (args_info.gce_given)
  {
    gce_ = true;
    gamma_.resize(args_info.gce_given);
    for (uint i=0; i!=args_info.gce_given; ++i)
      gamma_[i] = args_info.gce_arg[i];
  }

  if (args_info.structure_given)
  {
    data_str_list_.resize(args_info.structure_given);
    for (uint i=0; i!=args_info.structure_given; ++i)
      data_str_list_[i] = args_info.structure_arg[i];
  }

  if (args_info.reactivity_unpaired_given)
  {
    data_unpaired_list_.resize(args_info.reactivity_unpaired_given);
    for (uint i=0; i!=args_info.reactivity_unpaired_given; ++i)
      data_unpaired_list_[i] = args_info.reactivity_unpaired_arg[i];
  }

  if (args_info.reactivity_paired_given)
  {
    data_paired_list_.resize(args_info.reactivity_paired_given);
    for (uint i=0; i!=args_info.reactivity_paired_given; ++i)
      data_paired_list_[i] = args_info.reactivity_paired_arg[i];
  }

  noncomplementary_ = args_info.noncomplementary_flag==1;
  output_bpseq_ = args_info.bpseq_flag==1;
  pos_w_ = args_info.pos_w_arg;
  neg_w_ = args_info.neg_w_arg;
  lambda_ = args_info.lambda_arg;
  eta0_ = args_info.eta_arg;

  if ((!train_mode_ && args_info.inputs_num==0) ||
      (train_mode_ && data_str_list_.empty() && data_unpaired_list_.empty() && data_paired_list_.empty() && args_info.inputs_num==0)) 
  {
    cmdline_parser_print_help();
    cmdline_parser_free(&args_info);
    exit(1);
  }

  args_.resize(args_info.inputs_num);
  for (uint i=0; i!=args_info.inputs_num; ++i)
    args_[i]=args_info.inputs[i];

  cmdline_parser_free(&args_info);

  return *this;
}

std::pair<uint,uint>
read_data(std::vector<SStruct>& data, const std::vector<std::string>& lists, int type)
{
  std::pair<uint,uint> pos = std::make_pair(data.size(), data.size());
  for (auto l: lists)
  {
    std::ifstream is(l.c_str());
    if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + l);
    std::string f;
    while (is >> f)
      data.emplace_back(f, type);
  }
  pos.second = data.size();
  return pos;
}

int
NGSfold::train()
{
  uint n_first=1;

  // read traing data
  std::vector<SStruct> data;
  auto pos_str = read_data(data, data_str_list_, SStruct::NO_REACTIVITY);
  auto pos_unpaired = read_data(data, data_unpaired_list_, SStruct::REACTIVITY_UNPAIRED);
  auto pos_paired = read_data(data, data_paired_list_, SStruct::REACTIVITY_PAIRED);

  auto inference_engine = new InferenceEngine<double>(noncomplementary_);
  std::unique_ptr<ParameterHash<double>> pm(new ParameterHash<double>());
  if (!param_file_.empty())
    pm->ReadFromFile(param_file_);
  AdaGradRDAUpdater adagrad(eta0_, lambda_);

  // first round: train only from the full structure dataset
  std::vector<uint> idx(pos_str.second-pos_str.first);
  std::iota(idx.begin(), idx.end(), 0);
  for (uint epoch=0; epoch!=n_first; ++epoch)
  {
    std::random_shuffle(idx.begin(), idx.end());
    for (auto i : idx)
    {
      std::unordered_map<std::string,double> grad;
      inference_engine->LoadValues(std::move(pm));

      inference_engine->LoadSequence(data[i]);
      inference_engine->ComputeViterbi();
      auto pred = inference_engine->ComputeViterbiFeatureCounts();
      for (auto e : *pred)
        grad.insert(std::make_pair(e.first, 0.0)).first->second -= e.second;

      inference_engine->LoadSequence(data[i]);
      inference_engine->UseConstraints(data[i].GetMapping());
      inference_engine->ComputeViterbi();
      auto corr = inference_engine->ComputeViterbiFeatureCounts();
      for (auto e : *corr)
        grad.insert(std::make_pair(e.first, 0.0)).first->second += e.second;

      pm = inference_engine->LoadValues(nullptr);
      for (auto g : grad)
        adagrad.update(g.first, pm->get_by_key(g.first), g.second);
      adagrad.proceed_time();
    }
  }

  // second round

  delete inference_engine;

  return 0;
}

int
NGSfold::predict()
{
  // set parameters
  std::unique_ptr<ParameterHash<double>> pm(new ParameterHash<double>());

  if (!param_file_.empty())
    pm->ReadFromFile(param_file_);
  else if (noncomplementary_)
    pm->LoadFromHash(default_params_noncomplementary);
  else
    pm->LoadFromHash(default_params_complementary);
  
  // predict ss
  auto inference_engine = new InferenceEngine<double>(noncomplementary_);
  inference_engine->LoadValues(std::move(pm));
  for (auto s : args_)
  {
    SStruct sstruct;
    sstruct.Load(s);
    SStruct solution(sstruct);
    inference_engine->LoadSequence(sstruct);
    if (!mea_ && !gce_)
    {
      inference_engine->ComputeViterbi();
      solution.SetMapping(inference_engine->PredictPairingsViterbi());
    }
    else
    {
      inference_engine->ComputeInside();
      inference_engine->ComputeOutside();
      inference_engine->ComputePosterior();
      if (mea_)
        solution.SetMapping(inference_engine->PredictPairingsPosterior<0>(gamma_[0]));
      else
        solution.SetMapping(inference_engine->PredictPairingsPosterior<1>(gamma_[0]));
    }
    if (output_bpseq_)
      solution.WriteBPSEQ(std::cout);
    else
      solution.WriteParens(std::cout);
  }
  delete inference_engine;

  return 0;
}

int
main(int argc, char* argv[])
{
  try {
    return NGSfold().parse_options(argc, argv).run();
  } catch (const char* str) {
    std::cerr << str << std::endl;
  } catch (std::runtime_error e) {
    std::cerr << e.what() << std::endl;
  }
  return -1;
}
