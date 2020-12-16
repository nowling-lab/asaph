#!/usr/bin/env bats

setup() {
    N_INDIVIDUALS=20
    N_SNPS=10000
    
    export TEST_TEMP_DIR=`mktemp -u --tmpdir asaph-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}
    
    export VCF_PATH="${TEST_TEMP_DIR}/test.vcf"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    export PHENO_PATH="${TEST_TEMP_DIR}/phenotypes.txt"
    
    export CMD="asaph_generate_data"
}

@test "Run generate_data with no arguments" {
    run ${CMD}
    [ "$status" -eq 2 ]
}

@test "Run generate data with --help option" {
    run ${CMD} --help
}

@test "Test VCF data generation" {
    run ${CMD} \
        --n-populations 2 \
	--output-vcf ${VCF_PATH} \
	--output-populations ${POPS_PATH} \
	--individuals ${N_INDIVIDUALS} \
	--snps ${N_SNPS} \
        --n-phenotypes 3 \
        --output-phenotypes ${PHENO_PATH}

    [ "$status" -eq 0 ]
    [ -e ${VCF_PATH} ]
    [ -e ${POPS_PATH} ]
    [ -e ${PHENO_PATH} ]

    rm ${VCF_PATH}
    rm ${POPS_PATH}
    rm ${PHENO_PATH}
}

@test "Test VCF.gz data generation" {
    run ${CMD} \
        --n-populations 2 \
	--output-vcf-gz ${VCF_PATH}.gz \
	--output-populations ${POPS_PATH} \
	--individuals ${N_INDIVIDUALS} \
	--snps ${N_SNPS} \
        --n-phenotypes 3 \
        --output-phenotypes ${PHENO_PATH}

    [ "$status" -eq 0 ]
    [ -e ${VCF_PATH}.gz ]
    [ -e ${POPS_PATH} ]
    [ -e ${PHENO_PATH} ]

    rm ${VCF_PATH}.gz
    rm ${POPS_PATH}
    rm ${PHENO_PATH}
}
