#include <iostream>


#include "contrafold/Config.hpp"
#include "contrafold/Utilities.hpp"
#include "contrafold/InferenceEngine.hpp"
#include "contrafold/ParameterManager.hpp"
#include "contrafold/SStruct.hpp"
#include "contrafold/Defaults.ipp"

template 
class InferenceEngine<double>;

template
class ParameterManager<double>;

int main(int argc, const char* argv[])
{
  return 0;
}
