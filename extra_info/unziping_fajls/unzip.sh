#!/bin/bash

# Find all zip files in the current directory
zip_files=(*.zip)

# Loop through each zip file
for zip_file in "${zip_files[@]}"; do
    # Extract the zip file into a folder with the same name (without the .zip extension)
    folder_name=$(basename "$zip_file" .zip)
    mkdir -p "$folder_name" && unzip -d "$folder_name" "$zip_file"
done
