cmake_minimum_required(VERSION 3.16.3)

project(TiEMPO2 
    VERSION 0.0.2 
    DESCRIPTION "TiEMPO improved"
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(interputils SHARED src/include/InterpUtils.cpp)
target_include_directories(interputils PUBLIC src/include)

add_library(scan SHARED src/include/Scan.cpp)
target_include_directories(scan PUBLIC src/include)

add_library(tiempo2 SHARED src/cpp/InterfaceCPU.cpp)
target_link_libraries(tiempo2 PRIVATE interputils scan)
#target_include_directories(tiempo2 PRIVATE src/include)

if(NOT WIN32)
    target_compile_options(tiempo2 PRIVATE -pthread)
endif()

include(CheckLanguage)

check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 11)
    file(GLOB CUDAfiles src/cuda/*.cu)

    add_library(cutiempo2 SHARED ${CUDAfiles})
    target_include_directories(cutiempo2 PRIVATE src/cuda src/include)
    set_target_properties(cutiempo2 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET cutiempo2 PROPERTY CUDA_ARCHITECTURES OFF)
else()
    message(STATUS "WARNING: no CUDA compiler detected. Not building CUDA libraries")
endif()
