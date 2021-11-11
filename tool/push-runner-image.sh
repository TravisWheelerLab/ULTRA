#!/usr/bin/env sh

set -e

ULTRA_VERSION=$(./ultra -v)

docker push traviswheelerlab/ultra:${ULTRA_VERSION}
docker push traviswheelerlab/ultra:latest
