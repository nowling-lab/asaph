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

@test "Run snp_pair_associations workflow" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_pairwise_associations \
	--workdir ${WORKDIR_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations.txt" ]
}
