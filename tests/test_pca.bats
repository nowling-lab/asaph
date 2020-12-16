#!/usr/bin/env bats

load model_setup_helper

@test "Run pca with no arguments" {
    run asaph_pca
    [ "$status" -eq 2 ]
}

@test "Run pca with --help option" {
    run asaph_pca --help
    [ "$status" -eq 0 ]
}

@test "pca (categories)" {
    run asaph_pca \
        --workdir ${WORKDIR_PATH} \
        train \
        --n-components 6

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/pca.pkl" ]

    run asaph_pca \
	    --workdir ${WORKDIR_PATH} \
        explained-variance-analysis

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/pca_explained_variance_ratios.png" ]

    run asaph_pca \
	--workdir ${WORKDIR_PATH} \
        output-coordinates \
        --selected-components 1 2 3 4 \
        --output-fl ${WORKDIR_PATH}/pca_coordinates.txt

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.txt" ]
}

@test "pca (counts)" {
    run asaph_pca \
        --workdir ${COUNTS_WORKDIR_PATH} \
        train \
        --n-components 6

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/models/pca.pkl" ]

    run asaph_pca \
	    --workdir ${COUNTS_WORKDIR_PATH} \
        explained-variance-analysis

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/figures/pca_explained_variance_ratios.png" ]

    run asaph_pca \
	    --workdir ${COUNTS_WORKDIR_PATH} \
        output-coordinates \
        --selected-components 1 2 3 4 \
        --output-fl ${COUNTS_WORKDIR_PATH}/pca_coordinates.txt

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/pca_coordinates.txt" ]
}