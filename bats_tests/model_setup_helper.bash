setup() {
    N_INDIVIDUALS=20
    N_SNPS=250

    export TEST_TEMP_DIR=`mktemp -u --tmpdir asaph-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}

    export VCF_PATH="${TEST_TEMP_DIR}/test.vcf"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    export PHENO_PATH="${TEST_TEMP_DIR}/phenotypes.txt"
    export WORKDIR_PATH="${TEST_TEMP_DIR}/workdir"
    export COUNTS_WORKDIR_PATH="${TEST_TEMP_DIR}/counts_workdir"

    asaph_generate_data \
	--seed 1234 \
        --n-populations 2 \
	--output-vcf ${VCF_PATH} \
	--output-populations ${POPS_PATH} \
	--individuals ${N_INDIVIDUALS} \
	--snps ${N_SNPS} \
        --n-phenotypes 3 \
        --output-phenotypes ${PHENO_PATH}

    asaph_import \
	--workdir ${WORKDIR_PATH} \
	--vcf ${VCF_PATH} \
	--feature-type categories
    
    asaph_import \
	--workdir ${COUNTS_WORKDIR_PATH} \
	--vcf ${VCF_PATH} \
	--feature-type counts
}
