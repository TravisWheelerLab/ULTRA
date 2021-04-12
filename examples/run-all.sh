#!/usr/bin/env sh

set -e

run_one() {
  EX_PATH="$1"
  EX_LINE=$(echo "$EX_PATH" | tr -C "\n" "-")

  echo "+-$EX_LINE-+"
  echo "| $EX_PATH |"
  echo "+-$EX_LINE-+"
  echo ""
  ultra examples/ex1.fa
  echo ""
}

run_one examples/ex1.fa
run_one examples/ex2.fa
run_one examples/ex3.fa
run_one examples/ex4.fa
