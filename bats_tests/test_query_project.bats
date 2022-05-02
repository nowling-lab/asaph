#!/usr/bin/env bats

load model_setup_helper

file_contains() {
    local contains=`grep $2 $1 > /dev/null 2>&1 ; echo $?`
    echo $contains
}

@test "Run query_project with no arguments" {
    run asaph_query
    [ "$status" -eq 2 ]
}

@test "Run query_project with --help option" {
    run asaph_query --help
    [ "$status" -eq 0 ]
}

@test "Run query_project on test data" {
    RESULTS_FILE="${TEST_TEMP_DIR}/query_results.txt"
    run bash -c "asaph_query --workdir ${WORKDIR_PATH} > ${RESULTS_FILE}"

    [ "$status" -eq 0 ]
    [ $(file_contains ${RESULTS_FILE} ${N_INDIVIDUALS}) -eq 0 ]
}
