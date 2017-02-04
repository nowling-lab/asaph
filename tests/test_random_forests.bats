#!/usr/bin/env bats

@test "Run random_forests with no arguments" {
      
      run ${BATS_TEST_DIRNAME}/../bin/random_forests
      [ "$status" -eq 2 ]
}

@test "Run random_forests with --help option" {
      run ${BATS_TEST_DIRNAME}/../bin/random_forests --help
      [ "$status" -eq 0 ]
}