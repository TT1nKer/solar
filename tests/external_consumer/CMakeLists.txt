cmake_minimum_required(VERSION 3.20)

project(solar_external_consumer LANGUAGES CXX)

find_package(Solar CONFIG REQUIRED)

add_executable(solar_external_consumer probe.cpp)
target_compile_features(solar_external_consumer PRIVATE cxx_std_17)
target_link_libraries(solar_external_consumer PRIVATE Solar::Relativity)
