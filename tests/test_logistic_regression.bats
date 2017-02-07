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

@test "Train and output logistic regression model with sgd-l2 method" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--method sgd-l2

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr/sgd-l2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	rankings \
	--method sgd-l2

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings_lr_sgd-l2.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-l2.png" ]
}

@test "Train and output logistic regression model with sgd-en method" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
	--method sgd-en

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr/sgd-en" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	rankings \
	--method sgd-en

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings_lr_sgd-en.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-en.png" ]
}

@test "Train logistic regression model with default method" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr/sgd-l2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	rankings

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings_lr_sgd-l2.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-l2.png" ]
}
