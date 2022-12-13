
cmake_minimum_required(VERSION 3.12)
project(SSTreader C CXX)

find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED)

add_executable(SSTreader SSTreader.cpp)
target_link_libraries(SSTreader adios2::cxx11)

# if(ADIOS2_HAVE_MPI)
#   add_executable(hello_sstReader_mpi helloSstReader.cpp)
#   target_link_libraries(hello_sstReader_mpi adios2::cxx11_mpi MPI::MPI_C)
# endif()