#include <iostream>
#include <stdexcept>
#include "cmdline.h"
#include "Config.hpp"
#include "Utilities.hpp"
#include "InferenceEngine.hpp"
#include "ParameterHash.hpp"
#include "SStruct.hpp"

extern std::unordered_map<std::string, double> default_params_complementary;
extern std::unordered_map<std::string, double> default_params_noncomplementary;

class NGSfold
{
public:
  NGSfold() : train_mode_(false) { }

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

  if (args_info.train_given)
  {
    train_mode_ = true;
    out_file_ = args_info.train_arg;
  }
  else if (args_info.predict_given)
  {
    train_mode_ = false;
    param_file_ = args_info.predict_arg;
  }
  
  noncomplementary_ = args_info.noncomplementary_flag==1;
  output_bpseq_ = args_info.bpseq_flag==1;
  pos_w_ = args_info.pos_w_arg;
  neg_w_ = args_info.neg_w_arg;
  lambda_ = args_info.lambda_arg;
  eta0_ = args_info.eta_arg;

  if (args_info.inputs_num==0)
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

int
NGSfold::train()
{
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
  InferenceEngine<double>* inference_engine = new InferenceEngine<double>(noncomplementary_);
  SStruct sstruct;
  sstruct.Load(args_[0]);
  SStruct solution(sstruct);
  inference_engine->LoadValues(std::move(pm));
  inference_engine->LoadSequence(sstruct);
  inference_engine->ComputeViterbi();
  solution.SetMapping(inference_engine->PredictPairingsViterbi());
  delete inference_engine;

  if (output_bpseq_)
    solution.WriteBPSEQ(std::cout);
  else
    solution.WriteParens(std::cout);

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
