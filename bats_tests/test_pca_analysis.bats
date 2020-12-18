#!/usr/bin/env bats

load pca_setup_helper

@test "Run with no arguments" {
    run asaph_pca_analysis

    [ "$status" -eq 2 ]
}

@test "Run with --help option" {
    run asaph_pca_analysis --help

    [ "$status" -eq 0 ]
}

@test "plot projections (categories)" {
    run asaph_pca_analysis \
	--coordinates ${FULL_WORKDIR_PATH}/pca_coordinates.txt \
	plot-projections \
	--plot-dir ${FULL_WORKDIR_PATH}/pca_proj_plots \
	--pairs 1 2

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/pca_proj_plots/pca_projection_1_2.png" ]
}

@test "plot projections (hashed)" {
    run asaph_pca_analysis \
	--coordinates ${HASHED_WORKDIR_PATH}/pca_coordinates.txt \
	plot-projections \
	--plot-dir ${HASHED_WORKDIR_PATH}/pca_proj_plots \
	--pairs 1 2

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/pca_proj_plots/pca_projection_1_2.png" ]
}

@test "plot projections with labels (categories)" {
    run asaph_pca_analysis \
	--coordinates ${FULL_WORKDIR_PATH}/pca_coordinates.txt \
	plot-projections \
	--plot-dir ${FULL_WORKDIR_PATH}/pca_proj_plots \
	--pairs 1 2 \
	--labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/pca_proj_plots/pca_projection_1_2.png" ]
}

@test "plot projections (hashed)" {
    run asaph_pca_analysis \
	--coordinates ${HASHED_WORKDIR_PATH}/pca_coordinates.txt \
	plot-projections \
	--plot-dir ${HASHED_WORKDIR_PATH}/pca_proj_plots \
	--pairs 1 2 \
	--labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/pca_proj_plots/pca_projection_1_2.png" ]
}