project(hyperbolic-embedder)
cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 11)

add_executable(benchmark_embedder
        benchmark_embedder.cpp
        graph.cpp
        geometricds.cpp
        hyperbolic.cpp
        hyperboliclinear.cpp
        NLEHelper.cpp
        random.cpp
        springembedder.cpp
        plfit/powerlawCommon.cpp
        )

target_link_libraries(benchmark_embedder PRIVATE
        glog gflags gsl gslcblas)
