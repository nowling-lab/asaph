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

teardown() {
    rm -rf ${TEST_TEMP_DIR}
}
