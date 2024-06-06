#!/bin/bash

# Define experiment IDs
exp_ids=("deepmip-eocene-p1-PI" "deepmip-eocene-p1-x1" "deepmip-eocene-p1-x3" "deepmip-eocene-p1-x4")

# Define corresponding net files
net_files=(
    "COSMOS-landveg_r2413-piControl-rsnscs-v1.0.mean.nc"
    "COSMOS-landveg_r2413-deepmip_sens_1xCO2-rsnscs-v1.0.mean.nc"
    "COSMOS-landveg_r2413-deepmip_stand_3xCO2-rsnscs-v1.0.mean.nc"
    "COSMOS-landveg_r2413-deepmip_sens_4xCO2-rsnscs-v1.0.mean.nc"
)

# Loop over each experiment ID
for i in "${!exp_ids[@]}"; do
    exp="${exp_ids[$i]}"
    net="${net_files[$i]}"

    echo "Processing experiment: $exp"
    cd /Users/wb19586/Documents/coding_github/aprp_deepmip/data/deepmip-eocene-p1/COSMOS/COSMOS-landveg-r2413/$exp/v1.0/climatology

    down_all="rsds_COSMOS-landveg-r2413_${exp}_v1.0.mean.nc"
    down_all_masked="rsds_masked_COSMOS-landveg-r2413_${exp}_v1.0.mean.nc"
    down_mask="rsds_mask_COSMOS-landveg-r2413_${exp}_v1.0.mean.nc"

    up_all="rsus_COSMOS-landveg-r2413_${exp}_v1.0.mean.nc"

    down_cs="rsdscs_COSMOS-landveg-r2413_${exp}_v1.0.mean.nc"
    up_cs="rsuscs_COSMOS-landveg-r2413_${exp}_v1.0.mean.nc"

    albedo="albedo_${exp}.nc"

    echo "Processing experiment: $exp with net file: $net"

    # Set values to missing where down_all is between 0 and 0.1
    cdo setrtomiss,0,0.1 $down_all $down_all_masked

    # Calculate albedo
    cdo div $up_all $down_all_masked $albedo

    # Calculate down_cs
    cdo -chname,rsnscs,rsdscs -setmisstoc,0 -mulc,-1.0 -div $net -subc,1 $albedo $down_cs

    # Calculate up_cs
    cdo -chname,rsnscs,rsuscs -setmisstoc,0 -mulc,-1.0 -sub $net $down_cs $up_cs

    rm $down_all_masked $albedo

    echo "Completed processing for experiment: $exp"
done

echo "All experiments processed."