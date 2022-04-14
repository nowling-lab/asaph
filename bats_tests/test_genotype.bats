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

@test "plot projections (categories)" {
    run asaph_genotype \
	plot-projections \
	--workdir ${FULL_WORKDIR_PATH} \
	--pairs 1 2

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/plots/pca_projection_1_2.png" ]
}

@test "plot projections (hashed)" {
    run asaph_genotype \
	plot-projections \
	--workdir ${HASHED_WORKDIR_PATH} \
	--pairs 1 2

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/plots/pca_projection_1_2.png" ]
}

@test "plot projections with labels (categories)" {
    run asaph_genotype \
	plot-projections \
	--workdir ${FULL_WORKDIR_PATH} \
	--pairs 1 2 \
	--labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/plots/pca_projection_1_2.png" ]
}

@test "plot projections (hashed)" {
    run asaph_genotype \
	plot-projections \
	--workdir ${HASHED_WORKDIR_PATH} \
	--pairs 1 2 \
	--labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/plots/pca_projection_1_2.png" ]
}

@test "clustering (categories)" {
    run asaph_genotype \
    	unsupervised-genotyping \
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
    	unsupervised-genotyping \
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
