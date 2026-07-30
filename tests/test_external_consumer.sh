#!/usr/bin/env sh

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
solar_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
consumer_test_root=$(mktemp -d)

cleanup() {
    rm -rf -- "$consumer_test_root"
}
trap cleanup EXIT HUP INT TERM

solar_build_dir="$consumer_test_root/solar-build"
solar_install_dir="$consumer_test_root/solar-install"
consumer_build_dir="$consumer_test_root/consumer-build"

cmake \
    -S "$solar_root" \
    -B "$solar_build_dir" \
    -DSOLAR_BUILD_CLI=OFF \
    -DCMAKE_BUILD_TYPE=Release
cmake --build "$solar_build_dir" --parallel
cmake --install "$solar_build_dir" --prefix "$solar_install_dir"

cmake \
    -S "$solar_root/tests/external_consumer" \
    -B "$consumer_build_dir" \
    -DCMAKE_PREFIX_PATH="$solar_install_dir" \
    -DCMAKE_BUILD_TYPE=Release
cmake --build "$consumer_build_dir" --parallel
"$consumer_build_dir/solar_external_consumer"
