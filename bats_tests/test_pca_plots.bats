#!/usr/bin/env bats

load pca_setup_helper

@test "plot projections (categories)" {
    run asaph_pca \
	--workdir ${FULL_WORKDIR_PATH} \
	plot-projections \
	--pairs 1 2

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/plots/pca_projection_1_2.png" ]
}

@test "plot projections (hashed)" {
    run asaph_pca \
	--workdir ${HASHED_WORKDIR_PATH} \
	plot-projections \
	--pairs 1 2

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/plots/pca_projection_1_2.png" ]
}

@test "plot projections with labels (categories)" {
    run asaph_pca \
	--workdir ${FULL_WORKDIR_PATH} \
	plot-projections \
	--pairs 1 2 \
	--labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${FULL_WORKDIR_PATH}/plots/pca_projection_1_2.png" ]
}

@test "plot projections (hashed)" {
    run asaph_pca \
	--workdir ${HASHED_WORKDIR_PATH} \
	plot-projections \
	--pairs 1 2 \
	--labels-fl ${POPS_PATH}

    [ "$status" -eq 0 ]
    [ -e "${HASHED_WORKDIR_PATH}/plots/pca_projection_1_2.png" ]
}
