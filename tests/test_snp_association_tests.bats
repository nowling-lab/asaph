#!/usr/bin/env bats

load model_setup_helper

@test "Run snp_association_tests with no arguments" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests
    [ "$status" -eq 2 ]
}

@test "Run snp_association_tests with --help option" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests --help
    [ "$status" -eq 0 ]
}

@test "Calculate snp_association_tests (categories, default arguments)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, class probabilities intercept, adjusted training set)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --intercept none \
        --adjustment training-set

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, no intercept, adjusted training set)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --intercept none \
        --adjustment training-set

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, free-parameter intercept, adjusted training set)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --intercept free-parameter \
        --adjustment training-set

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, class probabilities intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --intercept class-probabilities \
        --adjustment none

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, no intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --intercept none \
        --adjustment none

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, free-parameter intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --intercept free-parameter \
        --adjustment none

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, remove empty columns)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --remove-empty-columns \
        --adjustment scaling-factor

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, adjusted scaling factor)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --adjustment scaling-factor

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}

@test "Calculate snp_association_tests (categories, pop dependent)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --populations ${POPS_PATH} \
        --dependent-variable genotype

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}
