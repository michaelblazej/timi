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

# Fetch Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(timi_tests test/test.cpp)
target_link_libraries(timi_tests PRIVATE timi GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(timi_tests)

# Examples
add_executable(basic_2d examples/basic_2d.cpp)
target_link_libraries(basic_2d PRIVATE timi)

add_executable(custom_class examples/custom_class.cpp)
target_link_libraries(custom_class PRIVATE timi)
