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
    --populations ${POPS_PATH} \
	--bagging \
	--method sgd-l2 \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method sgd-l2 \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
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
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-l2_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_sgd-l2_75.png" ]
}

@test "Logistic Regression workflow with no bagging, sgd-l2" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method sgd-l2 \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method sgd-l2 \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
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
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-l2_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_sgd-l2_75.png" ]
}

@test "Logistic Regression workflow with bagging, sgd-en" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method sgd-en \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method sgd-en \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
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
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-en_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-en_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_sgd-en_75.png" ]
}

@test "Logistic Regression workflow with no bagging, sgd-en" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method sgd-en \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method sgd-en \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-en/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
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
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-en_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-en_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_sgd-en_75.png" ]
}

@test "Logistic Regression workflow with bagging, sag-l2" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method sag-l2 \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method sag-l2 \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method sag-l2 \
	--n-models 150

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/150" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/150/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/150/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings \
	--method sag-l2

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sag-l2.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sag-l2.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--method sag-l2 \
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sag-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sag-l2_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_sag-l2_75.png" ]
}

@test "Logistic Regression workflow with no bagging, sag-l2" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method sag-l2 \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method sag-l2 \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method sag-l2 \
	--n-models 150

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/150" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/150/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sag-l2/150/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings \
	--method sag-l2

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sag-l2.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_sag-l2.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--method sag-l2 \
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sag-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sag-l2_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_sag-l2_75.png" ]
}

@test "Logistic Regression workflow with bagging, asgd-l2" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method asgd-l2 \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method asgd-l2 \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--method asgd-l2 \
	--n-models 150

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/150" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/150/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/150/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings \
	--method asgd-l2

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_asgd-l2.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_asgd-l2.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--method asgd-l2 \
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_asgd-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_asgd-l2_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_asgd-l2_75.png" ]
}

@test "Logistic Regression workflow with no bagging, asgd-l2" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method asgd-l2 \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method asgd-l2 \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--method asgd-l2 \
	--n-models 150

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/150" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/150/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-asgd-l2/150/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings \
	--method asgd-l2

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_asgd-l2.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_asgd-l2.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--method asgd-l2 \
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_asgd-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_asgd-l2_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_asgd-l2_75.png" ]
}

@test "Logistic Regression workflow with bagging, default method" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--bagging \
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
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-l2_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_sgd-l2_75.png" ]
}


@test "Logistic Regression workflow with no bagging, default method" {
    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--n-models 50

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/50/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
	--n-models 75

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/1" ]
    [ -e "${WORKDIR_PATH}/models/lr-sgd-l2/75/2" ]

    run ${BATS_TEST_DIRNAME}/../bin/logistic_regression \
	--workdir ${WORKDIR_PATH} \
	train \
    --populations ${POPS_PATH} \
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
	--n-models 75 \
	--plot-similarities

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings" ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_lr_sgd-l2_75.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/lr_weights_sgd-l2_75.png" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_lr_sgd-l2_75.png" ]
}
