cmake_minimum_required(VERSION 3.1)

#enable_language(cxx)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(MxFold)

find_package(PkgConfig REQUIRED)

pkg_check_modules(VIENNARNA REQUIRED RNAlib2)
include_directories(${VIENNARNA_INCLUDE_DIRS})
add_definitions(-DHAVE_VIENNA18)
add_definitions(-DHAVE_VIENNA20)

add_executable(
  mxfold
  src/main.cpp
  src/FeatureMap.cpp
  src/SStruct.cpp
  src/adagrad.cpp
  src/InferenceEngine.cpp
  src/Utilities.cpp
  src/default_params.cpp
  src/cmdline.c
  )
target_link_libraries(mxfold ${VIENNARNA_LDFLAGS})
