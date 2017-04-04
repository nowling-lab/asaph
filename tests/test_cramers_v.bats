#!/usr/bin/env bats

load model_setup_helper

@test "Run snp_associations with no arguments" {
    run ${BATS_TEST_DIRNAME}/../bin/cramers_v
    [ "$status" -eq 2 ]
}

@test "Run snp_associations with --help option" {
    run ${BATS_TEST_DIRNAME}/../bin/cramers_v --help
    [ "$status" -eq 0 ]
}

@test "Calculate Cramer's V vs population structure" {
    run ${BATS_TEST_DIRNAME}/../bin/cramers_v \
	    --workdir ${WORKDIR_PATH} \
        populations

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_population_associations.txt" ]
}

@test "Calculate Cramer's V pairwise" {
    run ${BATS_TEST_DIRNAME}/../bin/cramers_v \
	    --workdir ${WORKDIR_PATH} \
        pairwise

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations.txt" ]
}

@test "Calculate Cramer's V pairwise using sampling" {
    run ${BATS_TEST_DIRNAME}/../bin/cramers_v \
	    --workdir ${WORKDIR_PATH} \
        pairwise \
        --samples 100

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_pairwise_associations.txt" ]
}

@test "Calculate Cramer's V pairwise vs single SNP" {
    run ${BATS_TEST_DIRNAME}/../bin/cramers_v \
	    --workdir ${WORKDIR_PATH} \
        pairwise-single \
        --chrom 1 \
        --pos 200

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_associations_snp_1_200.txt" ]
}

