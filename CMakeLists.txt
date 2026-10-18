cmake_minimum_required(VERSION 3.16)
project(euler2d_mfem LANGUAGES CXX)

# Remove lint reports with the generated Makefile/Ninja clean target, including
# reports redirected into the build directory during out-of-source builds.
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/lint.log"
  "${CMAKE_CURRENT_BINARY_DIR}/lint.log"
)

find_package(MFEM CONFIG REQUIRED)

# Repair a stale dependency path present in some relocatable MFEM installations.
get_target_property(_mfem_links mfem INTERFACE_LINK_LIBRARIES)
if("/usr/lib64/libz.so" IN_LIST _mfem_links AND NOT EXISTS "/usr/lib64/libz.so")
  find_package(ZLIB REQUIRED)
  list(REMOVE_ITEM _mfem_links "/usr/lib64/libz.so")
  list(APPEND _mfem_links ZLIB::ZLIB)
  set_property(TARGET mfem PROPERTY INTERFACE_LINK_LIBRARIES "${_mfem_links}")
endif()

add_executable(euler2d
  src/main.cpp
  src/options.cpp
  src/physics.cpp
  src/problem.cpp
  src/operators.cpp
  src/partial.cpp
  src/condensation.cpp
  src/linear.cpp
  src/stepper.cpp
  src/verification.cpp
)
target_include_directories(euler2d PRIVATE include)
target_link_libraries(euler2d PRIVATE mfem)
target_compile_features(euler2d PRIVATE cxx_std_17)
if(MSVC)
  target_compile_options(euler2d PRIVATE /W4)
else()
  target_compile_options(euler2d PRIVATE -Wall -Wextra -Wpedantic)
endif()
# Opt in on the machine used for timings; preserve portable binaries by default.
option(EULER_NATIVE "Optimize for the build host CPU" OFF)
if(EULER_NATIVE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(euler2d PRIVATE -march=native)
endif()

# Optional existing algebra checks. Building does not execute them.
include(CTest)
if(BUILD_TESTING)
  add_test(NAME algebra COMMAND euler2d --verify)
endif()
