cmake_minimum_required(VERSION 3.0)

project(boost_geometry_simplify)

# set(CMAKE_CXX_FLAGS_DEBUG "-D_GLIBCXX_DEBUG")
FIND_PACKAGE( Boost 1.65 REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
ADD_EXECUTABLE(test test.cpp)
ADD_EXECUTABLE(example example.cpp)
ADD_EXECUTABLE(example_ls example_ls.cpp)
