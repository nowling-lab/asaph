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

@test "Calculate likelihood_ratio_test (categories, default arguments)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}

@test "Calculate likelihood_ratio_test (categories, class probabilities intercept, adjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH} \
        --intercept none \
        --training-set adjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}

@test "Calculate likelihood_ratio_test (categories, no intercept, adjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH} \
        --intercept none \
        --training-set adjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}

@test "Calculate likelihood_ratio_test (categories, free-parameter intercept, adjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH} \
        --intercept free-parameter \
        --training-set adjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}

@test "Calculate likelihood_ratio_test (categories, class probabilities intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH} \
        --intercept class-probabilities \
        --training-set unadjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}

@test "Calculate likelihood_ratio_test (categories, no intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH} \
        --intercept none \
        --training-set unadjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}

@test "Calculate likelihood_ratio_test (categories, free-parameter intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH} \
        --intercept free-parameter \
        --training-set unadjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}

@test "Calculate likelihood_ratio_test (categories, remove empty columns)" {
    run ${BATS_TEST_DIRNAME}/../bin/likelihood_ratio_test \
	    --workdir ${WORKDIR_PATH} \
        --remove-empty-columns

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_likelihood_ratio_tests.tsv" ]
}
