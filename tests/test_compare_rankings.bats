#!/usr/bin/env bats

@test "Run compare_rankings with no arguments" {
      
      run ${BATS_TEST_DIRNAME}/../bin/compare_rankings
      [ "$status" -eq 2 ]
}

@test "Run compare_rankings with --help option" {
      run ${BATS_TEST_DIRNAME}/../bin/compare_rankings --help
      [ "$status" -eq 0 ]
}