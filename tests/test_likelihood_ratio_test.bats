#!/usr/bin/env bats

load model_setup_helper

@test "Run likelihood_ratio_test with no arguments" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test
    [ "$status" -eq 2 ]
}

@test "Run likelihood_ratio_test with --help option" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test --help
    [ "$status" -eq 0 ]
}

@test "Calculate likelihood_ratio_test (categories)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}
