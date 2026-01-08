cmake_minimum_required(VERSION 3.14)
project(timi VERSION 1.0.0 LANGUAGES CXX)

# Header-only interface library
add_library(timi INTERFACE)
target_include_directories(timi INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(timi INTERFACE cxx_std_17)

# Tests
enable_testing()
add_executable(timi_tests test/test.cpp)
target_link_libraries(timi_tests PRIVATE timi)
add_test(NAME timi_tests COMMAND timi_tests)

# Examples
add_executable(basic_2d examples/basic_2d.cpp)
target_link_libraries(basic_2d PRIVATE timi)

add_executable(custom_class examples/custom_class.cpp)
target_link_libraries(custom_class PRIVATE timi)
