#!/usr/bin/env bats

load pca_setup_helper

@test "Run with no arguments" {
    run asaph_clustering

    [ "$status" -eq 2 ]
}

@test "Run with --help option" {
    run asaph_clustering --help

    [ "$status" -eq 0 ]
}

@test "kmeans clustering (categories)" {
    run asaph_clustering \
    	cluster-samples-kmeans \
	--coordinates ${FULL_WORKDIR_PATH}/pca_coordinates.txt \
	--n-clusters 3 \
	--components 1 2 \
	--output-labels-fl ${FULL_WORKDIR_PATH}/kmeans.labels

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/kmeans.labels" ]

    run asaph_clustering \
    	test-clusters \
	--cluster-labels-fl ${FULL_WORKDIR_PATH}/kmeans.labels \
	--other-labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
}

@test "dbscan clustering (categories)" {
    run asaph_clustering \
    	cluster-samples-dbscan \
	--coordinates ${FULL_WORKDIR_PATH}/pca_coordinates.txt \
	--components 1 \
	--output-labels-fl ${FULL_WORKDIR_PATH}/dbscan.labels

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/dbscan.labels" ]

    run asaph_clustering \
    	test-clusters \
	--cluster-labels-fl ${FULL_WORKDIR_PATH}/dbscan.labels \
	--other-labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
}

@test "kmeans clustering (hashed)" {
    run asaph_clustering \
    	cluster-samples-kmeans \
	--coordinates ${HASHED_WORKDIR_PATH}/pca_coordinates.txt \
	--n-clusters 3 \
	--components 1 2 \
	--output-labels-fl ${HASHED_WORKDIR_PATH}/kmeans.labels

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/kmeans.labels" ]

    run asaph_clustering \
    	test-clusters \
	--cluster-labels-fl ${HASHED_WORKDIR_PATH}/kmeans.labels \
	--other-labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
}

@test "dbscan clustering (hashed)" {
    run asaph_clustering \
    	cluster-samples-dbscan \
	--coordinates ${HASHED_WORKDIR_PATH}/pca_coordinates.txt \
	--components 1 \
	--epsilon 0.8 \
	--output-labels-fl ${HASHED_WORKDIR_PATH}/dbscan.labels

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/dbscan.labels" ]

    run asaph_clustering \
    	test-clusters \
	--cluster-labels-fl ${HASHED_WORKDIR_PATH}/dbscan.labels \
	--other-labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
}