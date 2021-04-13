cmake_minimum_required(VERSION 3.1.0)

file(GLOB ALLOTYPES_SRC *.cpp)

add_executable(AlloType ${ALLOTYPES_SRC})

target_link_libraries(AlloType Boost::filesystem Boost::system)
