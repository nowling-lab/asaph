#!/usr/bin/env bats

load import_helper

setup() {
    N_INDIVIDUALS=20
    N_SNPS=10000

    export TEST_TEMP_DIR=`mktemp -u --tmpdir asaph-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}

    export VCF_PATH="${TEST_TEMP_DIR}/test.vcf"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    export PHENO_PATH="${TEST_TEMP_DIR}/phenotypes.txt"
    export WORKDIR_PATH="${TEST_TEMP_DIR}/workdir"

    export IMPORT_CMD="asaph_pca"

    asaph_generate_data \
                        --seed 1234 \
                        --n-populations 2 \
                        --output-vcf ${VCF_PATH} \
                        --output-populations ${POPS_PATH} \
                        --individuals ${N_INDIVIDUALS} \
                        --snps ${N_SNPS} \
                        --n-phenotypes 3 \
                        --output-phenotypes ${PHENO_PATH}

    gzip -k ${VCF_PATH}
}

@test "Run import with no arguments" {
    run ${IMPORT_CMD}
    [ "$status" -eq 2 ]
}

@test "Run import with --help option" {
    run ${IMPORT_CMD} --help
    [ "$status" -eq 0 ]
}

@test "Import data: vcf, default" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, counts" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --feature-type allele-counts

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, categories" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --feature-type genotype-categories

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, feature hashing" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --feature-type allele-counts \
	    --sampling-method feature-hashing

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, counts" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type allele-counts

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, categories" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type genotype-categories

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, feature hashing, reduced by size" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type allele-counts \
	    --sampling-method feature-hashing \
	    --min-inversion-fraction 0.05

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, feature hashing, reduced by dimensions" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type allele-counts \
	    --sampling-method feature-hashing \
	    --num-dimensions 10

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_features ${WORKDIR_PATH}) -eq 10 ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, categories, reservoir, reduced by size" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type genotype-categories \
	    --sampling-method reservoir \
	    --min-inversion-fraction 0.05

    N_FEATURES=$((N_SNPS * 3))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, categories, reservoir, reduced by dimensions" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type genotype-categories \
	    --sampling-method reservoir \
	    --num-dimensions 10

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ -e "${WORKDIR_PATH}/pca_coordinates.tsv" ]
    [ $(count_features ${WORKDIR_PATH}) -eq 10 ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}
