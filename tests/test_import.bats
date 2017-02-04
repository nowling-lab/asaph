#!/usr/bin/env bats

@test "Run import with no arguments" {
      
      run ${BATS_TEST_DIRNAME}/../bin/import
      [ "$status" -eq 2 ]
}

@test "Run import with --help option" {
      run ${BATS_TEST_DIRNAME}/../bin/import --help
      [ "$status" -eq 0 ]
}