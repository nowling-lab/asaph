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
	    --workdir ${WORKDIR_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}

@test "Calculate snp_association_tests (categories, class probabilities intercept, adjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --intercept none \
        --training-set adjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}

@test "Calculate snp_association_tests (categories, no intercept, adjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --intercept none \
        --training-set adjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}

@test "Calculate snp_association_tests (categories, free-parameter intercept, adjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --intercept free-parameter \
        --training-set adjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}

@test "Calculate snp_association_tests (categories, class probabilities intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --intercept class-probabilities \
        --training-set unadjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}

@test "Calculate snp_association_tests (categories, no intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --intercept none \
        --training-set unadjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}

@test "Calculate snp_association_tests (categories, free-parameter intercept, unadjusted)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --intercept free-parameter \
        --training-set unadjusted

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}

@test "Calculate snp_association_tests (categories, remove empty columns)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --remove-empty-columns

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_gt.tsv" ]
}

@test "Calculate snp_association_tests (categories, pop dependent)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --dependent-variable population

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_lrtests_pop.tsv" ]
}
