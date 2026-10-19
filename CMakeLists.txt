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
  src/block_ilu.cpp
  src/stepper.cpp
  src/verification.cpp
)
target_include_directories(euler2d PRIVATE include)
target_link_libraries(euler2d PRIVATE mfem)

# MFEM build-tree package exports created with MFEM_FETCH_TPLS currently name
# HYPRE/METIS in their link interface but do not export the fetched include and
# library directories.  Installed MFEM packages do not need this workaround;
# make the in-tree configuration used by build-intel-mfem.sh usable as well.
get_target_property(_mfem_includes mfem INTERFACE_INCLUDE_DIRECTORIES)
foreach(_mfem_include IN LISTS _mfem_includes)
  if(EXISTS "${_mfem_include}/fetch/hypre/include/HYPRE_utilities.h")
    target_include_directories(euler2d SYSTEM PRIVATE
      "${_mfem_include}/fetch/hypre/include"
      "${_mfem_include}/fetch/metis/include")
    target_link_directories(euler2d PRIVATE
      "${_mfem_include}/fetch/hypre/lib"
      "${_mfem_include}/fetch/metis/lib")
    if(EXISTS "${_mfem_include}/fetch/metis/lib/libmetis.a")
      get_target_property(_mfem_build_links mfem INTERFACE_LINK_LIBRARIES)
      list(REMOVE_ITEM _mfem_build_links METIS)
      list(APPEND _mfem_build_links "${_mfem_include}/fetch/metis/lib/libmetis.a")
      set_property(TARGET mfem PROPERTY INTERFACE_LINK_LIBRARIES
        "${_mfem_build_links}")
    endif()
  endif()
endforeach()
target_compile_features(euler2d PRIVATE cxx_std_17)

# Hot residual/operator kernels cross several translation units. Enable supported
# Release-mode interprocedural optimization so these small calls can be inlined
# without changing floating-point semantics (unlike fast-math).
include(CheckIPOSupported)
check_ipo_supported(RESULT EULER_IPO_SUPPORTED OUTPUT EULER_IPO_ERROR)
if(EULER_IPO_SUPPORTED)
  set_property(TARGET euler2d PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

if(MSVC)
  target_compile_options(euler2d PRIVATE /W4)
else()
  target_compile_options(euler2d PRIVATE -Wall -Wextra -Wpedantic)
endif()
# Opt in on the machine used for timings; preserve portable binaries by default.
option(EULER_NATIVE "Optimize for the build host CPU" OFF)
if(EULER_NATIVE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang|Intel|IntelLLVM)$")
    target_compile_options(euler2d PRIVATE -march=native)
  endif()
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  target_compile_options(euler2d PRIVATE -qopenmp-simd)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(euler2d PRIVATE -fopenmp-simd)
endif()

# Optional existing algebra checks. Building does not execute them.
include(CTest)
if(BUILD_TESTING)
  add_test(NAME algebra COMMAND euler2d --verify)
endif()
