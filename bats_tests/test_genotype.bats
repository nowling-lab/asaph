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
	--coordinates ${FULL_WORKDIR_PATH}/pca_coordinates.txt \
	--plot-dir ${FULL_WORKDIR_PATH}/pca_proj_plots \
	--pairs 1 2

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/pca_proj_plots/pca_projection_1_2.png" ]
}

@test "plot projections (hashed)" {
    run asaph_genotype \
	plot-projections \
	--coordinates ${HASHED_WORKDIR_PATH}/pca_coordinates.txt \
	--plot-dir ${HASHED_WORKDIR_PATH}/pca_proj_plots \
	--pairs 1 2

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/pca_proj_plots/pca_projection_1_2.png" ]
}

@test "plot projections with labels (categories)" {
    run asaph_genotype \
	plot-projections \
	--coordinates ${FULL_WORKDIR_PATH}/pca_coordinates.txt \
	--plot-dir ${FULL_WORKDIR_PATH}/pca_proj_plots \
	--pairs 1 2 \
	--labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/pca_proj_plots/pca_projection_1_2.png" ]
}

@test "plot projections (hashed)" {
    run asaph_genotype \
	plot-projections \
	--coordinates ${HASHED_WORKDIR_PATH}/pca_coordinates.txt \
	--plot-dir ${HASHED_WORKDIR_PATH}/pca_proj_plots \
	--pairs 1 2 \
	--labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/pca_proj_plots/pca_projection_1_2.png" ]
}

@test "dbscan clustering (categories)" {
    run asaph_genotype \
    	cluster-samples \
	--coordinates ${FULL_WORKDIR_PATH}/pca_coordinates.txt \
	--components 1 \
	--output-labels-fl ${FULL_WORKDIR_PATH}/dbscan.labels

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/dbscan.labels" ]

    run asaph_genotype \
    	evaluate-clusters \
	--cluster-labels-fl ${FULL_WORKDIR_PATH}/dbscan.labels \
	--other-labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
}

@test "dbscan clustering (hashed)" {
    run asaph_genotype \
    	cluster-samples \
	--coordinates ${HASHED_WORKDIR_PATH}/pca_coordinates.txt \
	--components 1 \
	--epsilon 0.8 \
	--output-labels-fl ${HASHED_WORKDIR_PATH}/dbscan.labels

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/dbscan.labels" ]

    run asaph_genotype \
    	evaluate-clusters \
	--cluster-labels-fl ${HASHED_WORKDIR_PATH}/dbscan.labels \
	--other-labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
}