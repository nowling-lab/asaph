#!/usr/bin/env bats

load pca_setup_helper

@test "Run with no arguments" {
    run asaph_genotype

    [ "$status" -eq 2 ]
}

@test "Run with --help option" {
    run asaph_genotype --help

    [ "$status" -eq 0 ]
}

@test "clustering (categories)" {
    run asaph_genotype \
    	cluster \
	--workdir ${FULL_WORKDIR_PATH} \
	--components 1 \
	--n-clusters 3 \
	--predicted-labels-fl ${FULL_WORKDIR_PATH}/unsupervised.labels

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/unsupervised.labels" ]

    run asaph_genotype \
    	evaluate-predicted-genotypes \
	--predicted-labels-fl ${FULL_WORKDIR_PATH}/unsupervised.labels \
	--known-labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
}

@test "clustering (hashed)" {
    run asaph_genotype \
    	cluster \
	--workdir ${HASHED_WORKDIR_PATH} \
	--components 1 \
	--n-clusters 3 \
	--predicted-labels-fl ${HASHED_WORKDIR_PATH}/unsupervised.labels

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/unsupervised.labels" ]

    run asaph_genotype \
    	evaluate-predicted-genotypes \
	--predicted-labels-fl ${HASHED_WORKDIR_PATH}/unsupervised.labels \
	--known-labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
}
