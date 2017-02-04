#!/usr/bin/env bats


@test "Run logistic_regression with no arguments" {
      
      run ${BATS_TEST_DIRNAME}/../bin/logistic_regression
      [ "$status" -eq 2 ]
}

@test "Run logistic_regression with --help option" {
      run ${BATS_TEST_DIRNAME}/../bin/logistic_regression --help
      [ "$status" -eq 0 ]
}