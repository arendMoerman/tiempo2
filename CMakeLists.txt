cmake_minimum_required(VERSION 3.16.3)

project(TiEMPO2 
    VERSION 0.0.1 
    DESCRIPTION "TiEMPO improved"
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(interputils SHARED src/include/InterpUtils.cpp)
target_include_directories(interputils PUBLIC src/include)

add_library(scan SHARED src/include/Scan.cpp)
target_include_directories(scan PUBLIC src/include)

add_library(observe SHARED src/include/Observe.cpp)
target_include_directories(observe PUBLIC src/include)

add_library(tiempo2 SHARED src/cpp/InterfaceCPU.cpp)
target_link_libraries(tiempo2 PRIVATE interputils scan observe)
#target_include_directories(tiempo2 PRIVATE src/include)

if(NOT WIN32)
    target_compile_options(tiempo2 PRIVATE -pthread)
endif()

