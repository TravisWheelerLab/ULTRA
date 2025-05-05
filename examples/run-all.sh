#!/usr/bin/env sh

set -e

run_one() {
  EX_PATH="$1"
  EX_LINE=$(echo "$EX_PATH" | tr -C "\n" "-")

  echo "+-$EX_LINE-+"
  echo "| $EX_PATH |"
  echo "+-$EX_LINE-+"
  echo ""
  ./ultra "$EX_PATH"
  echo ""
}

run_one examples/example_1.fa
run_one examples/example_2.fa
run_one examples/example_3.fa