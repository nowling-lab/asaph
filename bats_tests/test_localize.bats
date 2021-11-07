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
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]

    run asaph_detect_and_localize \
        association-tests \
        --pca-coordinates-tsv ${WORKDIR_PATH}/pca_coordinates.tsv \
	    --vcf ${VCF_PATH} \
	    --components 1 2 \
	    --pca-associations-tsv ${TEST_TEMP_DIR}/pca_tests.tsv

    [ "$status" -eq 0 ]
    [ -e "${TEST_TEMP_DIR}/pca_tests.tsv" ]

    run asaph_detect_and_localize \
        plot \
        --pca-associations-tsv "${TEST_TEMP_DIR}/pca_tests.tsv" \
	    --plot-fl "${TEST_TEMP_DIR}/manhattan_plot_comp1.png" \
	    --component 1 \
	    --chromosome 1 \
	    --n-windows -1

   [ "$status" -eq 0 ]
   [ -e "${TEST_TEMP_DIR}/manhattan_plot_comp1.png" ]

    run asaph_detect_and_localize \
        plot \
        --pca-associations-tsv "${TEST_TEMP_DIR}/pca_tests.tsv" \
        --plot-fl "${TEST_TEMP_DIR}/manhattan_plot_comp2.png" \
	    --component 2 \
	    --n-windows -1

   [ "$status" -eq 0 ]
   [ -e "${TEST_TEMP_DIR}/manhattan_plot_comp2.png" ]
}

@test "pca (counts)" {
    [ -e "${COUNTS_WORKDIR_PATH}/models/model" ]
    [ -e "${COUNTS_WORKDIR_PATH}/pca_coordinates.tsv" ]

    run asaph_detect_and_localize \
        association-tests \
        --pca-coordinates-tsv ${COUNTS_WORKDIR_PATH}/pca_coordinates.tsv \
	    --vcf ${VCF_PATH} \
	    --components 1 2 \
	    --pca-associations-tsv ${TEST_TEMP_DIR}/pca_tests_counts.tsv

    [ "$status" -eq 0 ]
    [ -e "${TEST_TEMP_DIR}/pca_tests_counts.tsv" ]

    run asaph_detect_and_localize \
        plot \
        --pca-associations-tsv "${TEST_TEMP_DIR}/pca_tests_counts.tsv" \
	    --plot-fl "${TEST_TEMP_DIR}/manhattan_plot_comp1.png" \
	    --component 1 \
	    --n-windows -1

   [ "$status" -eq 0 ]
   [ -e "${TEST_TEMP_DIR}/manhattan_plot_comp1.png" ]

    run asaph_detect_and_localize \
        plot \
        --pca-associations-tsv "${TEST_TEMP_DIR}/pca_tests_counts.tsv" \
	    --plot-fl "${TEST_TEMP_DIR}/manhattan_plot_comp2.png" \
	    --component 2 \
	    --n-windows -1

   [ "$status" -eq 0 ]
   [ -e "${TEST_TEMP_DIR}/manhattan_plot_comp2.png" ]
}
