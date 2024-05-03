#!/bin/bash

# Directory containing the .nc files
SOURCE_DIR="../aprp_output_data/deepmip"

# Subdirectory to store the remapped files
TARGET_DIR="${SOURCE_DIR}/remapped"

# Create the target directory if it doesn't exist
mkdir -p "${TARGET_DIR}"

# Define the grid file
GRID="r180x90"  # Update this with the path to your grid file

# Loop through all .nc files in the source directory
for FILE in ${SOURCE_DIR}/*.nc; do
    # Extract the filename without the path and extension
    BASENAME=$(basename -- "$FILE")
    FILENAME="${BASENAME%.nc}"

    # Define the output filename
    OUTPUT_FILE="${TARGET_DIR}/${FILENAME}.r180x90.nc"

    # Remap the file
    cdo sethalo,0,1 -remapbil,${GRID} "${FILE}" "${OUTPUT_FILE}"
done

echo "Remapping completed."
