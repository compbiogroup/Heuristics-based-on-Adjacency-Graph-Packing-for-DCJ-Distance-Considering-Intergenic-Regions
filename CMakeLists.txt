cmake_minimum_required(VERSION 3.5)

project(dec_ir)

option(CXX "enable C++ compilation" ON)
enable_language(CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

file(GLOB MISC misc/*.hpp misc/*.cpp)
file(GLOB EXTER external/*.hpp external/*.cpp)
file(GLOB CYCLE cycle/*.hpp cycle/*.cpp)
file(GLOB HEUR heur/*.hpp heur/*.cpp)

set(CMAKE_CXX_STANDARD 14)
set(CXX_FILESYSTEM_LIBRARIES "stdc++fs")

#################################################################################################
# decomposition
#################################################################################################

find_package(OpenMP REQUIRED)
add_executable(dec main.cpp ${MISC} ${CYCLE} ${HEUR} ${GRIMM})
target_compile_options(dec PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(dec PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

#################################################################################################

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()
