#!/usr/bin/env sh

set -e

cmake .
make
ULTRA_VERSION=$(./ultra -v)

docker build -f Dockerfile_run \
    -t traviswheelerlab/ultra:${ULTRA_VERSION} \
    -t traviswheelerlab/ultra:latest \
    $@ \
    .
