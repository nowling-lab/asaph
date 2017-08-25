setup() {
    N_INDIVIDUALS=20
    N_SNPS=1000

    export TEST_TEMP_DIR=`mktemp -u --tmpdir asaph-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}

    export VCF_PATH="${TEST_TEMP_DIR}/test.vcf"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    export WORKDIR_PATH="${TEST_TEMP_DIR}/workdir"
    export COUNTS_WORKDIR_PATH="${TEST_TEMP_DIR}/counts_workdir"

    ${BATS_TEST_DIRNAME}/../bin/generate_data \
			            --seed 1234 \
                        --n-populations 2 \
			            --output-vcf ${VCF_PATH} \
			            --output-populations ${POPS_PATH} \
			            --individuals ${N_INDIVIDUALS} \
			            --snps ${N_SNPS}

    ${BATS_TEST_DIRNAME}/../bin/import \
			            --workdir ${WORKDIR_PATH} \
			            --vcf ${VCF_PATH} \
			            --populations ${POPS_PATH} \
			            --feature-type categories
    
    ${BATS_TEST_DIRNAME}/../bin/import \
			            --workdir ${COUNTS_WORKDIR_PATH} \
			            --vcf ${VCF_PATH} \
			            --populations ${POPS_PATH} \
			            --feature-type counts
}

teardown() {
    rm -rf ${TEST_TEMP_DIR}
}
