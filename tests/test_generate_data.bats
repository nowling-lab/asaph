#!/usr/bin/env bats

setup() {
    N_INDIVIDUALS=20
    N_SNPS=10000
    
    export TEST_TEMP_DIR=`dirname $(mktemp -u)`
    
    export VCF_PATH="${TEST_TEMP_DIR}/test.vcf"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    
    export CMD="${BATS_TEST_DIRNAME}/../bin/generate_data"
}

@test "Run generate_data with no arguments" {
    run ${CMD}
    [ "$status" -eq 2 ]
}

@test "Run generate data with --help option" {
    run ${CMD} --help
}

@test "Test data generation" {
    run ${CMD} \
	--output-vcf ${VCF_PATH} \
	--output-populations ${POPS_PATH} \
	--individuals ${N_INDIVIDUALS} \
	--snps ${N_SNPS}

    [ "$status" -eq 0 ]
    [ -e ${VCF_PATH} ]
    [ -e ${POPS_PATH} ]

    rm ${VCF_PATH}
    rm ${POPS_PATH}
}
