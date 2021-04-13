cmake_minimum_required(VERSION 3.1.0)

project(AlloType)

set(CMAKE_CXX_STANDARD   11)
set(CMAKE_CXX_STANDARD_REQUIRED   ON)

# Find Boost library
set(Boost_USE_MULTITHREADED     ON)
set(Boost_USE_STATIC_LIBS       ON)
set(Boost_USE_STATIC_RUNTIME   OFF)

find_package(Boost COMPONENTS filesystem system REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIR})
endif()

# Find Eigen3 library
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(src)

install(TARGETS AlloType RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
