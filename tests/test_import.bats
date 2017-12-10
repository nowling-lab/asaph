#!/usr/bin/env bats

load import_helper

setup() {
    N_INDIVIDUALS=20
    N_SNPS=10000

    export TEST_TEMP_DIR=`mktemp -u --tmpdir asaph-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}

    export VCF_PATH="${TEST_TEMP_DIR}/test.vcf"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    export WORKDIR_PATH="${TEST_TEMP_DIR}/workdir"

    export IMPORT_CMD="${BATS_TEST_DIRNAME}/../bin/import"

    ${BATS_TEST_DIRNAME}/../bin/generate_data \
                        --seed 1234 \
                        --n-populations 2 \
                        --output-vcf ${VCF_PATH} \
                        --output-populations ${POPS_PATH} \
                        --individuals ${N_INDIVIDUALS} \
                        --snps ${N_SNPS}
}

@test "Run import with no arguments" {
    run ${IMPORT_CMD}
    [ "$status" -eq 2 ]
}

@test "Run import with --help option" {
    run ${IMPORT_CMD} --help
    [ "$status" -eq 0 ]
}

@test "Import data: vcf, dna, counts" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --populations ${POPS_PATH} \
	    --feature-type counts

    N_FEATURE_INDICES=$((N_SNPS))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/snp_feature_genotypes" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
    [ $(count_snp_feature_indices ${WORKDIR_PATH}) -eq ${N_FEATURE_INDICES} ]
}

@test "Import data: vcf, dna, categories" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --populations ${POPS_PATH} \
	    --feature-type categories

    N_FEATURE_INDICES=$((N_SNPS * 3))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/snp_feature_genotypes" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
    [ $(count_snp_feature_indices ${WORKDIR_PATH}) -eq ${N_FEATURE_INDICES} ]
}

@test "Import data: vcf, dna, counts, feature_compression" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --populations ${POPS_PATH} \
	    --feature-type counts \
	    --compress

    N_FEATURE_INDICES=$((N_SNPS))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
    [ $(count_snp_feature_indices ${WORKDIR_PATH}) -eq ${N_FEATURE_INDICES} ]
}

@test "Import data: vcf, dna, categories, feature_compression" {
    run ${IMPORT_CMD} \
	    --workdir ${WORKDIR_PATH} \
	    --vcf ${VCF_PATH} \
	    --populations ${POPS_PATH} \
	    --feature-type categories \
	    --compress

    N_FEATURE_INDICES=$((N_SNPS * 3))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ -e "${WORKDIR_PATH}/snp_feature_genotypes" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
    [ $(count_snp_feature_indices ${WORKDIR_PATH}) -eq ${N_FEATURE_INDICES} ]
}
