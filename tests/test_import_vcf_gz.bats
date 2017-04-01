#!/usr/bin/env bats

load import_helper

@test "Import data: vcf.gz, dna, counts" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf-gz ${VCF_GZ_PATH} \
	--populations ${POPS_PATH} \
	--feature-type counts

    N_FEATURE_INDICES=$((N_SNPS * 2))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
    [ $(count_snp_feature_indices ${WORKDIR_PATH}) -eq ${N_FEATURE_INDICES} ]
}

@test "Import data: vcf.gz, dna, categories" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf-gz ${VCF_GZ_PATH} \
	--populations ${POPS_PATH} \
	--feature-type categories

    N_FEATURE_INDICES=$((N_SNPS * 3))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
    [ $(count_snp_feature_indices ${WORKDIR_PATH}) -eq ${N_FEATURE_INDICES} ]
}

@test "Import data: vcf.gz, dna, counts, feature_compression" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf-gz ${VCF_GZ_PATH} \
	--populations ${POPS_PATH} \
	--feature-type counts \
	--compress

    N_FEATURE_INDICES=$((N_SNPS * 2))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
    [ $(count_snp_feature_indices ${WORKDIR_PATH}) -eq ${N_FEATURE_INDICES} ]
}

@test "Import data: vcf.gz, dna, categories, feature_compression" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf-gz ${VCF_GZ_PATH} \
	--populations ${POPS_PATH} \
	--feature-type categories \
	--compress

    N_FEATURE_INDICES=$((N_SNPS * 3))

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ -e "${WORKDIR_PATH}/project_summary" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
    [ $(count_snp_feature_indices ${WORKDIR_PATH}) -eq ${N_FEATURE_INDICES} ]
}
