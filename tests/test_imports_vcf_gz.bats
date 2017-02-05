#!/usr/bin/env bats

count_snps() {
    local counts=`python -c "import cPickle; data=cPickle.load(open('${1}/snp_feature_indices')); print len(data)"`
    echo "$counts"
}

setup() {
    N_INDIVIDUALS=20
    N_SNPS=10000

    export TEST_TEMP_DIR=`dirname $(mktemp -u)`
    
    export VCF_GZ_PATH="${TEST_TEMP_DIR}/test.vcf.gz"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    export WORKDIR_PATH="${TEST_TEMP_DIR}/workdir"
    
    export IMPORT_CMD="${BATS_TEST_DIRNAME}/../bin/import"
    
    ${BATS_TEST_DIRNAME}/../bin/generate_data \
			--output-vcf-gz ${VCF_GZ_PATH} \
			--output-populations ${POPS_PATH} \
			--individuals ${N_INDIVIDUALS} \
			--snps ${N_SNPS}
}

teardown() {
    rm ${VCF_GZ_PATH}
    rm ${POPS_PATH}
    rm -rf ${WORKDIR_PATH}
}

@test "Import data: vcf.gz, dna, counts" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf-gz ${VCF_GZ_PATH} \
	--populations ${POPS_PATH} \
	--feature-type counts
    
    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
}

@test "Import data: vcf.gz, dna, categories" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf-gz ${VCF_GZ_PATH} \
	--populations ${POPS_PATH} \
	--feature-type categories

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
}

@test "Import data: vcf.gz, dna, counts, feature_compression" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf-gz ${VCF_GZ_PATH} \
	--populations ${POPS_PATH} \
	--feature-type counts \
	--compress

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
}

@test "Import data: vcf.gz, dna, categories, feature_compression" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf-gz ${VCF_GZ_PATH} \
	--populations ${POPS_PATH} \
	--feature-type categories \
	--compress

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
}
