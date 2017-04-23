count_snps() {
    local counts=`python -c "import cPickle; data=cPickle.load(open('${1}/snp_feature_indices')); print len(data)"`
    echo "$counts"
}

count_samples() {
    local counts=`python -c "import cPickle; data=cPickle.load(open('${1}/sample_labels')); print len(data)"`
    echo "$counts"
}

count_snp_feature_indices() {
    local count=`python -c "import cPickle; data=cPickle.load(open('${1}/snp_feature_indices')); print sum([len(idx) for idx in data.itervalues()])"`
    echo "$count"
}

setup() {
    N_INDIVIDUALS=20
    N_SNPS=10000

    export TEST_TEMP_DIR=`mktemp -u --tmpdir asaph-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}

    export VCF_GZ_PATH="${TEST_TEMP_DIR}/test.vcf.gz"
    export POPS_PATH="${TEST_TEMP_DIR}/populations.txt"
    export WORKDIR_PATH="${TEST_TEMP_DIR}/workdir"

    export IMPORT_CMD="${BATS_TEST_DIRNAME}/../bin/import"

    ${BATS_TEST_DIRNAME}/../bin/generate_data \
                        --seed 1234 \
                        --output-vcf-gz ${VCF_GZ_PATH} \
                        --output-populations ${POPS_PATH} \
                        --individuals ${N_INDIVIDUALS} \
                        --snps ${N_SNPS}
}

teardown() {
    rm -rf ${TEST_TEMP_DIR}
}
