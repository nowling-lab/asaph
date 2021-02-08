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
    [ -e "${WORKDIR_PATH}/models/model" ]

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

    run asaph_pca_association_tests \
        --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --components 1 2 \
	    --output-tsv ${TEST_TEMP_DIR}/pca_tests.tsv

    [ "$status" -eq 0 ]
    [ -e "${TEST_TEMP_DIR}/pca_tests.tsv" ]

    run asaph_detect_and_localize \
        plot \
        --input-tsv "${TEST_TEMP_DIR}/pca_tests.tsv" \
	    --plot-fl "${TEST_TEMP_DIR}/manhattan_plot_comp1.png" \
	    --component 1 \
	    --chromosome 1 \
	    --n-windows -1

   [ "$status" -eq 0 ]
   [ -e "${TEST_TEMP_DIR}/manhattan_plot_comp1.png" ]

    run asaph_detect_and_localize \
        plot \
        --input-tsv "${TEST_TEMP_DIR}/pca_tests.tsv" \
        --plot-fl "${TEST_TEMP_DIR}/manhattan_plot_comp2.png" \
	    --component 2 \
	    --chromosome 1 \
	    --n-windows -1

   [ "$status" -eq 0 ]
   [ -e "${TEST_TEMP_DIR}/manhattan_plot_comp2.png" ]
}

@test "pca (counts)" {
    [ -e "${COUNTS_WORKDIR_PATH}/models/model" ]

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

    run asaph_pca_association_tests \
        --workdir ${COUNTS_WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --components 1 2 \
	    --output-tsv ${TEST_TEMP_DIR}/pca_tests_counts.tsv

    [ "$status" -eq 0 ]
    [ -e "${TEST_TEMP_DIR}/pca_tests_counts.tsv" ]

    run asaph_detect_and_localize \
        plot \
        --input-tsv "${TEST_TEMP_DIR}/pca_tests_counts.tsv" \
	    --plot-fl "${TEST_TEMP_DIR}/manhattan_plot_comp1.png" \
	    --component 1 \
	    --chromosome 1 \
	    --n-windows -1

   [ "$status" -eq 0 ]
   [ -e "${TEST_TEMP_DIR}/manhattan_plot_comp1.png" ]

    run asaph_detect_and_localize \
        plot \
        --input-tsv "${TEST_TEMP_DIR}/pca_tests_counts.tsv" \
	    --plot-fl "${TEST_TEMP_DIR}/manhattan_plot_comp2.png" \
	    --component 2 \
        --chromosome 1 \
	    --n-windows -1

   [ "$status" -eq 0 ]
   [ -e "${TEST_TEMP_DIR}/manhattan_plot_comp2.png" ]
}
