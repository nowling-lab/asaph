#!/usr/bin/env bats

count_snps() {
    local counts=`python -c "import cPickle; data=cPickle.load(open('${1}/snp_feature_indices')); print len(data)"`
    echo "$counts"
}

count_samples() {
    local counts=`python -c "import cPickle; data=cPickle.load(open('${1}/sample_labels')); print len(data)"`
    echo "$counts"
}

setup() {
    N_INDIVIDUALS=20
    N_SNPS=10000

    export TEST_TEMP_DIR=`dirname $(mktemp -u)`

    export VCF_PATH="${TEST_TEMP_DIR}/test.vcf"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    export WORKDIR_PATH="${TEST_TEMP_DIR}/workdir"

    export IMPORT_CMD="${BATS_TEST_DIRNAME}/../bin/import"

    ${BATS_TEST_DIRNAME}/../bin/generate_data \
			--seed 1234 \
			--output-vcf ${VCF_PATH} \
			--output-populations ${POPS_PATH} \
			--individuals ${N_INDIVIDUALS} \
			--snps ${N_SNPS}
}

teardown() {
    rm ${VCF_PATH}
    rm ${POPS_PATH}
    rm -rf ${WORKDIR_PATH}
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

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, dna, categories" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf ${VCF_PATH} \
	--populations ${POPS_PATH} \
	--feature-type categories

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, dna, counts, feature_compression" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf ${VCF_PATH} \
	--populations ${POPS_PATH} \
	--feature-type counts \
	--compress

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}

@test "Import data: vcf, dna, categories, feature_compression" {
    run ${IMPORT_CMD} \
	--workdir ${WORKDIR_PATH} \
	--vcf ${VCF_PATH} \
	--populations ${POPS_PATH} \
	--feature-type categories \
	--compress

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}" ]
    [ -d "${WORKDIR_PATH}" ]
    [ $(count_snps ${WORKDIR_PATH}) -eq ${N_SNPS} ]
    [ $(count_samples ${WORKDIR_PATH}) -eq ${N_INDIVIDUALS} ]
}
