#!/usr/bin/env bats

load model_setup_helper

@test "Run pca with no arguments" {
    run ${BATS_TEST_DIRNAME}/../bin/pca
    [ "$status" -eq 2 ]
}

@test "Run pca with --help option" {
    run ${BATS_TEST_DIRNAME}/../bin/pca --help
    [ "$status" -eq 0 ]
}

@test "pca (categories)" {
    run ${BATS_TEST_DIRNAME}/../bin/pca \
        --workdir ${WORKDIR_PATH} \
        train \
        --n-components 6

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/pca.pkl" ]
    
    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${WORKDIR_PATH} \
        explained-variance-analysis

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/pca_explained_variance_ratios.png" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${WORKDIR_PATH} \
        plot-projections \
        --pairs 0 1 2 3

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/pca_projection_0_1.png" ]
    [ -e "${WORKDIR_PATH}/figures/pca_projection_2_3.png" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${WORKDIR_PATH} \
        output-coordinates \
        --selected-components 0 1 2 3 \
        --output-fl ${WORKDIR_PATH}/pca_coordinates.txt

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.txt" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${WORKDIR_PATH} \
        analyze-weights \
        --weights 1.0 1.0

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/pca_feature_weights_pc0.png" ]
    [ -e "${WORKDIR_PATH}/figures/pca_feature_weights_pc1.png" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${WORKDIR_PATH} \
        extract-genotypes \
        --weights 1.0 1.0 \
        --thresholds 0.5 0.5

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/analysis/component_genotypes.tsv" ]
}

@test "pca (counts)" {
    run ${BATS_TEST_DIRNAME}/../bin/pca \
        --workdir ${COUNTS_WORKDIR_PATH} \
        train \
        --n-components 6

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/models/pca.pkl" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${COUNTS_WORKDIR_PATH} \
        explained-variance-analysis

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/figures/pca_explained_variance_ratios.png" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${COUNTS_WORKDIR_PATH} \
        plot-projections \
        --pairs 0 1 2 3

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/figures/pca_projection_0_1.png" ]
    [ -e "${COUNTS_WORKDIR_PATH}/figures/pca_projection_2_3.png" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${COUNTS_WORKDIR_PATH} \
        output-coordinates \
        --selected-components 0 1 2 3 \
        --output-fl ${COUNTS_WORKDIR_PATH}/pca_coordinates.txt

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/pca_coordinates.txt" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${COUNTS_WORKDIR_PATH} \
        analyze-weights \
        --weights 1.0 1.0

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/figures/pca_feature_weights_pc0.png" ]
    [ -e "${COUNTS_WORKDIR_PATH}/figures/pca_feature_weights_pc1.png" ]

    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${COUNTS_WORKDIR_PATH} \
        extract-genotypes \
        --weights 1.0 1.0 \
        --thresholds 0.5 0.5

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/analysis/component_genotypes.tsv" ]
}

@test "Find min components to achieve explained variance threshold (categories)" {
    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${WORKDIR_PATH} \
        min-components-explained-variance \
        --init-n-components 10 \
        --explained-variance-threshold 0.05

    [ "$status" -eq 0 ]
}

@test "Find min components to achieve explained variance threshold (counts)" {
    run ${BATS_TEST_DIRNAME}/../bin/pca \
	    --workdir ${COUNTS_WORKDIR_PATH} \
        min-components-explained-variance \
        --init-n-components 10 \
        --explained-variance-threshold 0.05

    [ "$status" -eq 0 ]
}

