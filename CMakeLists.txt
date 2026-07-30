cmake_minimum_required(VERSION 3.20)

project(
    Solar
    VERSION 0.2.0
    DESCRIPTION "Solar-system dynamics and relativistic ray-tracing foundations"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(SOLAR_BUILD_CLI "Build the Solar command-line application" ON)

file(GLOB SOLAR_CORE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB SOLAR_RELATIVITY_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/relativity/*.cpp")

add_library(solar_core ${SOLAR_CORE_SOURCES})
add_library(Solar::Core ALIAS solar_core)
set_target_properties(solar_core PROPERTIES EXPORT_NAME Core)
target_compile_features(solar_core PUBLIC cxx_std_17)
target_include_directories(
    solar_core
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

add_library(solar_relativity ${SOLAR_RELATIVITY_SOURCES})
add_library(Solar::Relativity ALIAS solar_relativity)
set_target_properties(solar_relativity PROPERTIES EXPORT_NAME Relativity)
target_compile_features(solar_relativity PUBLIC cxx_std_17)
target_link_libraries(solar_relativity PUBLIC Solar::Core)

if(SOLAR_BUILD_CLI)
    add_executable(
        solar_cli
        cli/main.cpp
        cli/relativity_metric.cpp)
    set_target_properties(solar_cli PROPERTIES OUTPUT_NAME solar)
    target_compile_features(solar_cli PRIVATE cxx_std_17)
    target_link_libraries(solar_cli PRIVATE Solar::Relativity)
endif()

install(
    TARGETS solar_core solar_relativity
    EXPORT SolarTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

if(SOLAR_BUILD_CLI)
    install(TARGETS solar_cli RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()

install(
    DIRECTORY include/solar
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

set(SOLAR_CMAKE_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/Solar")

configure_package_config_file(
    cmake/SolarConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/SolarConfig.cmake"
    INSTALL_DESTINATION "${SOLAR_CMAKE_INSTALL_DIR}")

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SolarConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion)

install(
    EXPORT SolarTargets
    FILE SolarTargets.cmake
    NAMESPACE Solar::
    DESTINATION "${SOLAR_CMAKE_INSTALL_DIR}")
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/SolarConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/SolarConfigVersion.cmake"
    DESTINATION "${SOLAR_CMAKE_INSTALL_DIR}")
