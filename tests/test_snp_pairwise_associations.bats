#!/usr/bin/env bats

load model_setup_helper

@test "Run snp_pair_associations with no arguments" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations
    [ "$status" -eq 2 ]
}

@test "Run snp_pair_associations with --help option" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations --help
    [ "$status" -eq 0 ]
}

@test "Run snp_pair_associations workflow with defaults" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations \
	    --workdir ${WORKDIR_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations.txt" ]
}

@test "Run snp_pair_associations workflow with permutation model" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations \
	    --workdir ${WORKDIR_PATH} \
        --model-type permutation

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations_permutation.txt" ]
}

@test "Run snp_pair_associations workflow with uniform-random model" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations \
	    --workdir ${WORKDIR_PATH} \
        --model-type uniform-random

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations_uniform-random.txt" ]
}

@test "Run snp_pair_associations workflow with defaults, sampling" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations \
	    --workdir ${WORKDIR_PATH} \
        --samples 100

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations.txt" ]
}

@test "Run snp_pair_associations workflow with samples & permutation model" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations \
	    --workdir ${WORKDIR_PATH} \
        --model-type permutation \
        --samples 100

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations_permutation.txt" ]
}

@test "Run snp_pair_associations workflow with samples & uniform-random model" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations \
	    --workdir ${WORKDIR_PATH} \
        --model-type uniform-random \
        --samples 100

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations_uniform-random.txt" ]
}
