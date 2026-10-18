#!/bin/sh
# Static source checks only: no compilation, solver execution, or benchmark.
set -eu
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
clang-format --dry-run --Werror "$project_dir"/include/*.hpp "$project_dir"/src/*.cpp
cppcheck --enable=warning,style,performance,portability --std=c++17 \
  --inline-suppr --suppress=missingIncludeSystem --error-exitcode=1 \
  -I "$project_dir/include" "$project_dir/src"
