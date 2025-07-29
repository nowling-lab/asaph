#!/usr/bin/env bats

load model_setup_helper

@test "Run with no arguments" {
    run asaph_localize
    [ "$status" -eq 2 ]
}

@test "Run with --help option" {
    run asaph_localize --help
    [ "$status" -eq 0 ]
}

@test "pca (categories)" {
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]

    run asaph_localize \
        --workdir ${WORKDIR_PATH} \
        association-tests \
	--vcf ${VCF_PATH} \
	--components 1 2

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/pca_associations.tsv" ]

    run asaph_localize \
	--workdir ${WORKDIR_PATH} \
        manhattan-plot \
	--component 1 \
	--chromosome 1 \
	--y-limit 10

   [ "$status" -eq 0 ]
   [ -e "${WORKDIR_PATH}/plots/manhattan_pc1_chrom1.png" ]

    run asaph_localize \
        --workdir "${WORKDIR_PATH}" \
	manhattan-plot \
	--component 2 \
	--chromosome 1 \
	--y-limit 10

   [ "$status" -eq 0 ]
   [ -e "${WORKDIR_PATH}/plots/manhattan_pc2_chrom1.png" ]
}

@test "pca (counts)" {
    [ -e "${COUNTS_WORKDIR_PATH}/models/model" ]
    [ -e "${COUNTS_WORKDIR_PATH}/pca_coordinates.tsv" ]

    run asaph_localize \
	--workdir "${COUNTS_WORKDIR_PATH}" \
        association-tests \
	--vcf ${VCF_PATH} \
	--components 1 2

    [ "$status" -eq 0 ]
    [ -e "${COUNTS_WORKDIR_PATH}/pca_associations.tsv" ]

    run asaph_localize \
        --workdir ${COUNTS_WORKDIR_PATH} \
	manhattan-plot \
	--component 1 \
	--chromosome 1 \
	--y-limit 10

   [ "$status" -eq 0 ]
   [ -e "${COUNTS_WORKDIR_PATH}/plots/manhattan_pc1_chrom1.png" ]

    run asaph_localize \
        --workdir "${COUNTS_WORKDIR_PATH}" \
	manhattan-plot \
	--chromosome 1 \
	--component 2 \
	--y-limit 10

   [ "$status" -eq 0 ]
   [ -e "${COUNTS_WORKDIR_PATH}/plots/manhattan_pc2_chrom1.png" ]
}
