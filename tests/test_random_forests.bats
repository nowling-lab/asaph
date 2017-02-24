#!/usr/bin/env bats

load model_setup_helper

@test "Run random_forests with no arguments" {
    run ${BATS_TEST_DIRNAME}/../bin/random_forests
    [ "$status" -eq 2 ]
}

@test "Run random_forests with --help option" {
    run ${BATS_TEST_DIRNAME}/../bin/random_forests --help
    [ "$status" -eq 0 ]
}

@test "Run random_forests workflow" {
    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	train \
	--statistics \
	--interactions \
	--trees 100

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/rf/100" ]
    [ -d "${WORKDIR_PATH}/models/rf/100" ]
    [ -e "${WORKDIR_PATH}/models/rf/100/interactions" ]
    [ -e "${WORKDIR_PATH}/figures/features_used_histogram_rf_100_trees.png" ]

    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	train \
	--trees 250

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/rf/250" ]
    [ -d "${WORKDIR_PATH}/models/rf/250" ]

    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	train \
	--trees 500

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/rf/500" ]
    [ -d "${WORKDIR_PATH}/models/rf/500" ]

    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_rf.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_rf.pdf" ]
    [ -e "${WORKDIR_PATH}/figures/snp_counts.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_counts.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--trees 100

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_rf_100.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_rf_100.png" ]
}

@test "Run random_forests workflow with resampling" {
    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	train \
	--trees 100 \
	--statistics \
	--interactions \
	--resamples 10

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/rf/100" ]
    [ -d "${WORKDIR_PATH}/models/rf/100" ]
    [ -e "${WORKDIR_PATH}/models/rf/100/interactions" ]
    [ -e "${WORKDIR_PATH}/figures/features_used_histogram_rf_100_trees.png" ]

    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	train \
	--trees 250 \
	--resamples 10

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/rf/250" ]
    [ -d "${WORKDIR_PATH}/models/rf/250" ]

    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	train \
	--trees 500 \
	--resamples 10

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/models/rf/500" ]
    [ -d "${WORKDIR_PATH}/models/rf/500" ]

    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	analyze-rankings

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_rf.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_ranking_overlaps_rf.pdf" ]
    [ -e "${WORKDIR_PATH}/figures/snp_counts.png" ]
    [ -e "${WORKDIR_PATH}/figures/snp_counts.pdf" ]

    run ${BATS_TEST_DIRNAME}/../bin/random_forests \
	--workdir ${WORKDIR_PATH} \
	output-rankings \
	--trees 100

    [ "$status" -eq 0 ]
    [ -e "${WORKDIR_PATH}/rankings/rankings_rf_100.tsv" ]
    [ -e "${WORKDIR_PATH}/figures/within_model_similarity_rf_100.png" ]
}
