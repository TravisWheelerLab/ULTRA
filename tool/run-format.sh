#!/usr/bin/env sh

set -e

CPP_FILES=$(find src -type f -name '*.cpp')
HPP_FILES=$(find src -type f -name '*.hpp')

clang-format -i ${CPP_FILES} ${HPP_FILES}
