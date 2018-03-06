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

@test "Calculate snp_association_tests (categories)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_association_tests.tsv" ]
}

@test "Calculate snp_association_tests (categories, variables file)" {
    run ${BATS_TEST_DIRNAME}/../bin/snp_association_tests \
	    --workdir ${WORKDIR_PATH} \
        --variables-fl ${PHENO_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/statistics/snp_association_tests.tsv" ]
}
