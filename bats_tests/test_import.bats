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
	    train \
	    --vcf ${VCF_PATH}

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, counts" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf ${VCF_PATH} \
	    --feature-type counts

    N_FEATURES=$((N_SNPS * 2))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_features ${WORKDIR_PATH}) -eq ${N_FEATURES} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, categories" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf ${VCF_PATH} \
	    --feature-type categories

    N_FEATURES=$((N_SNPS * 3))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_features ${WORKDIR_PATH}) -eq ${N_FEATURES} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, hashed" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf ${VCF_PATH} \
	    --feature-type hashed

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, counts" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type counts

    N_FEATURES=$((N_SNPS * 2))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_features ${WORKDIR_PATH}) -eq ${N_FEATURES} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, categories" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type categories

    N_FEATURES=$((N_SNPS * 3))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_features ${WORKDIR_PATH}) -eq ${N_FEATURES} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, hashed, reduced by size" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type hashed \
	    --min-inversion-fraction 0.05

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, hashed, reduced by dimensions" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type hashed \
	    --num-dimensions 10

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_features ${WORKDIR_PATH}) -eq 10 ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, categories, reservoir, reduced by size" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type categories \
	    --subsampling-method reservoir \
	    --min-inversion-fraction 0.05

    N_FEATURES=$((N_SNPS * 3))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, categories, reservoir, reduced by dimensions" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type categories \
	    --subsampling-method reservoir \
	    --num-dimensions 10

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_features ${WORKDIR_PATH}) -eq 10 ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, hashed, random projection, reduced by size" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type hashed \
	    --subsampling-method random-projection \
	    --min-inversion-fraction 0.2 \
	    --inner-dim 4096

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf.gz, hashed, random-projection, reduced by dimensions" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    train \
	    --vcf-gz ${VCF_PATH}.gz \
	    --feature-type hashed \
	    --subsampling-method random-projection \
	    --num-dimensions 10 \
	    --inner-dim 1024

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/models/model" ]
    [ $(count_features ${WORKDIR_PATH}) -eq 10 ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}