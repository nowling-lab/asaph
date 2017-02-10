#!/usr/bin/env bats

load model_setup_helper

@test "Run logistic_regression with no arguments" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression
    [ "$status" -eq 2 ]
}

@test "Run logistic_regression with --help option" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression --help
    [ "$status" -eq 0 ]
}

@test "Logistic Regression workflow with bagging, sgd-l2" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--method sgd-l2 \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--method sgd-l2 \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--method sgd-l2 \
	--n-models 150

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/150" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/150/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/150/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings \
	--method sgd-l2

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sgd-l2.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sgd-l2.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--method sgd-l2 \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-l2_75.png" ]
}

@test "Logistic Regression workflow with bagging, sgd-en" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--method sgd-en \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--method sgd-en \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--method sgd-en \
	--n-models 150

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/150" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/150/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/150/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings \
	--method sgd-en

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sgd-en.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sgd-en.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--method sgd-en \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-en_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-en_75.png" ]
}

@test "Logistic Regression workflow with bagging, default method" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/2" ]

        run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--n-models 150

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/150" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/150/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/150/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sgd-l2.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sgd-l2.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-l2_75.png" ]
}
