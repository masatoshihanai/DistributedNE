cmake_minimum_required(VERSION 3.5)
project(DistributedNE)

set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -D_GLIBCXX_PARALLEL -fopenmp")

add_compile_options(-g)
add_compile_options(-DNDEBUG)

###### User defined parameters #####
## UINT32_VERTT: Vertex ID type becomes 32-bit unsigned int instead of 64-bit
#add_compile_options(-DUINT32_VERTT)
#
## MAX_NUM_PARTS: The maximum number of partitions. Default is 1024
#add_compile_options(-DMAX_NUM_PARTS=2028)
#
####################################

find_library(MPI Names mpi PATHS)
message(STATUS "${emph}mpi ${emphoff} ${MPI}")

include_directories(. include ${CMAKE_BINARY_DIR}/include)
add_executable(DistributedNE src/main.cpp)
target_link_libraries(DistributedNE ${MPI})
