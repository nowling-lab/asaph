count_features() {
    local counts=`asaph_query --workdir ${1} | grep n_features | cut -d ' ' -f 2`
    echo "$counts"
}

count_samples() {
    local counts=`asaph_query --workdir ${1} | grep n_samples | cut -d ' ' -f 2`
    echo "$counts"
}
