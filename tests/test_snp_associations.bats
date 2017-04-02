#!/usr/bin/env bats

load model_setup_helper

@test "Run snp_associations with no arguments" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_associations
    [ "$status" -eq 2 ]
}

@test "Run snp_associations with --help option" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_associations --help
    [ "$status" -eq 0 ]
}

@test "Run snp_associations workflow with SNP" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_associations \
	    --workdir ${WORKDIR_PATH} \
        snp \
        --chrom 1 \
        --pos 300

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_associations_snp_1_300.txt" ]
}
