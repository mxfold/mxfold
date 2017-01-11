#include "../config.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <string>
#include <stdexcept>
#include <ctime>
#include "../config.h"
#include "cmdline.h"
#include "Config.hpp"
#include "Utilities.hpp"
#include "InferenceEngine.hpp"
#include "ParameterHash.hpp"
#include "SStruct.hpp"
#include "adagrad.hpp"

extern std::unordered_map<std::string, double> default_params_complementary;
extern std::unordered_map<std::string, double> default_params_noncomplementary;

class NGSfold
{
public:
  NGSfold() : train_mode_(false), mea_(false), gce_(false) { }

  NGSfold& parse_options(int& argc, char**& argv);

  int run()
  {
    if (validation_mode_)
      return validate();
    if (train_mode_)
      return train();
    else
      return predict();
  }

private:
  int train();
  int predict();
  int validate();
  std::pair<uint,uint> read_data(std::vector<SStruct>& data, const std::vector<std::string>& lists, int type) const;
  std::unordered_map<std::string,double> compute_gradients(const SStruct& s, InferenceEngine<double>* inference_engine);
  

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
  uint t_max_;
  uint t_burn_in_;
  float weight_weak_labeled_;
  float pos_w_;
  float neg_w_;
  float lambda_;
  float eta0_;
  float scale_reactivity_;
  float threshold_unpaired_reactivity_;
  float threshold_paired_reactivity_;
  bool discretize_reactivity_;
  int verbose_;
  std::string out_param_;
  bool validation_mode_;
  bool use_constraints_;
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

  if (args_info.unpaired_reactivity_given)
  {
    data_unpaired_list_.resize(args_info.unpaired_reactivity_given);
    for (uint i=0; i!=args_info.unpaired_reactivity_given; ++i)
      data_unpaired_list_[i] = args_info.unpaired_reactivity_arg[i];
  }

  if (args_info.paired_reactivity_given)
  {
    data_paired_list_.resize(args_info.paired_reactivity_given);
    for (uint i=0; i!=args_info.paired_reactivity_given; ++i)
      data_paired_list_[i] = args_info.paired_reactivity_arg[i];
  }

  if (args_info.out_param_given)
    out_param_ = args_info.out_param_arg;

  noncomplementary_ = args_info.noncomplementary_flag==1;
  output_bpseq_ = args_info.bpseq_flag==1;
  t_max_ = args_info.max_iter_arg;
  t_burn_in_ = args_info.burn_in_arg;
  weight_weak_labeled_ = args_info.weight_weak_label_arg;
  pos_w_ = args_info.pos_w_arg;
  neg_w_ = args_info.neg_w_arg;
  lambda_ = args_info.lambda_arg;
  eta0_ = args_info.eta_arg;
  scale_reactivity_ = args_info.scale_reactivity_arg;
  threshold_unpaired_reactivity_ = args_info.threshold_unpaired_reactivity_arg;
  threshold_paired_reactivity_ = args_info.threshold_paired_reactivity_arg;
  discretize_reactivity_ = args_info.discretize_reactivity_flag==1;
  verbose_ = args_info.verbose_arg;
  use_constraints_ = args_info.constraints_flag==1;
  validation_mode_ = args_info.validate_flag==1;

  srand(args_info.random_seed_arg<0 ? time(0) : args_info.random_seed_arg);

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
NGSfold::
read_data(std::vector<SStruct>& data, const std::vector<std::string>& lists, int type) const
{
  std::pair<uint,uint> pos = std::make_pair(data.size(), data.size());
  for (auto l: lists)
  {
    std::ifstream is(l.c_str());
    if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + l);
    std::string f;
    while (is >> f)
    {
      data.emplace_back(f, type);
      if (type!=SStruct::NO_REACTIVITY && discretize_reactivity_)
        data.back().DiscretizeReactivity(threshold_unpaired_reactivity_, threshold_paired_reactivity_);
    }
  }
  pos.second = data.size();
  return pos;
}

std::unordered_map<std::string,double>
NGSfold::
compute_gradients(const SStruct& s, InferenceEngine<double>* inference_engine)
{
  double starting_time = GetSystemTime();
  std::unordered_map<std::string,double> grad;
  // count the occurence of parameters in the predicted structure
  inference_engine->LoadSequence(s);
  if (s.GetType() == SStruct::NO_REACTIVITY)
    inference_engine->UseLossBasePair(s.GetMapping(), pos_w_, neg_w_);
  else if (discretize_reactivity_)
    inference_engine->UseLossPosition(s.GetMapping(), pos_w_, neg_w_);
  else
    inference_engine->UseLossReactivity(s.GetReactivityUnpair(), s.GetReactivityPair(), pos_w_, neg_w_);
  inference_engine->ComputeViterbi();
  auto loss1 = inference_engine->GetViterbiScore();
  if (verbose_>1)
  {
    SStruct solution(s);
    solution.SetMapping(inference_engine->PredictPairingsViterbi());
    std::cout << std::endl;
    solution.WriteParens(std::cout);
  }
  auto pred = inference_engine->ComputeViterbiFeatureCounts();
  for (auto e : *pred)
    grad.insert(std::make_pair(e.first, 0.0)).first->second += e.second;

  // count the occurence of parameters in the correct structure
  if (s.GetType() == SStruct::NO_REACTIVITY || discretize_reactivity_)
  {
    inference_engine->LoadSequence(s);
    inference_engine->UseConstraints(s.GetMapping());
  }
  else
  {
    inference_engine->LoadSequence(s, true, threshold_unpaired_reactivity_, threshold_paired_reactivity_, scale_reactivity_);
  }    
  inference_engine->ComputeViterbi();
  auto loss2 = inference_engine->GetViterbiScore();
  if (verbose_>1)
  {
    SStruct solution(s);
    solution.SetMapping(inference_engine->PredictPairingsViterbi());
    solution.WriteParens(std::cout);
  }
  auto corr = inference_engine->ComputeViterbiFeatureCounts();
  for (auto e : *corr)
    grad.insert(std::make_pair(e.first, 0.0)).first->second -= e.second;

  if (verbose_>0)
  {
    std::cout << "Seq: " << s.GetNames()[0] << ", "
              << "Loss: " << loss1-loss2 << ", "
              << "Time: " << GetSystemTime() - starting_time << "sec" << std::endl;
  }

  return std::move(grad);
}

int
NGSfold::train()
{
  // read traing data
  std::vector<SStruct> data;
  auto pos_str = read_data(data, data_str_list_, SStruct::NO_REACTIVITY);
  auto pos_unpaired = read_data(data, data_unpaired_list_, SStruct::REACTIVITY_UNPAIRED);
  auto pos_paired = read_data(data, data_paired_list_, SStruct::REACTIVITY_PAIRED);

  // set up the inference engine
  auto inference_engine = new InferenceEngine<double>(noncomplementary_);
  std::unique_ptr<ParameterHash<double>> pm(new ParameterHash<double>());
  if (!param_file_.empty())
    pm->ReadFromFile(param_file_);
  AdaGradRDAUpdater adagrad(eta0_, lambda_);

  // run max-margin training
  for (uint t=0, k=0; t!=t_max_; ++t)
  {
    if (verbose_>0)
      std::cout << std::endl << "=== Epoch " << t << " ===" << std::endl;
    std::vector<uint> idx(t<t_burn_in_ ? pos_str.second : pos_paired.second);
    std::iota(idx.begin(), idx.end(), 0);
    std::random_shuffle(idx.begin(), idx.end());
    for (auto i : idx)
    {
      auto w = i<pos_str.second ? 1.0 : weight_weak_labeled_;
      inference_engine->LoadValues(std::move(pm));
      auto grad = compute_gradients(data[i], inference_engine);
      pm = inference_engine->LoadValues(nullptr);
      for (auto g : grad)
        if (g.second!=0.0)
          adagrad.update(g.first, pm->get_by_key(g.first), g.second*w);
      adagrad.proceed_time();

      if (!out_param_.empty())
        pm->WriteToFile(SPrintF("%s/%d.param", out_param_.c_str(), k++));
    }
  }

  delete inference_engine;
  
  pm->WriteToFile(out_file_);

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
    //sstruct.WriteParens(std::cout); // for debug
    //std::cout << std::endl;
    if (use_constraints_)
      inference_engine->UseConstraints(sstruct.GetMapping());
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
NGSfold::validate()
{
  // set parameters
  std::unique_ptr<ParameterHash<double>> pm(new ParameterHash<double>());

  if (!param_file_.empty())
    pm->ReadFromFile(param_file_);
  else if (noncomplementary_)
    pm->LoadFromHash(default_params_noncomplementary);
  else
    pm->LoadFromHash(default_params_complementary);
  
  for (auto s : args_)
  {
    SStruct sstruct;
    sstruct.Load(s);
    SStruct solution(sstruct);
    InferenceEngine<double> inference_engine(noncomplementary_, 
                                             std::max<int>(sstruct.GetLength()/5., DEFAULT_C_MAX_SINGLE_LENGTH));
    inference_engine.LoadValues(std::move(pm));
    inference_engine.LoadSequence(sstruct);
    inference_engine.UseConstraints(sstruct.GetMapping());
    inference_engine.ComputeViterbi();
    auto score = inference_engine.GetViterbiScore();
    std::cout << sstruct.GetNames()[0] << " " 
              << (score>NEG_INF ? "OK" : "NG") << std::endl;
    if (score==NEG_INF)
    {
      sstruct.WriteParens(std::cout); // for debug
      std::cout << std::endl;
    }  
    pm = inference_engine.LoadValues(nullptr);
  }

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
