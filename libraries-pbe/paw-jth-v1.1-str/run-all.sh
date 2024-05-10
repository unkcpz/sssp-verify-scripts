#!/bin/bash

input_folder="."

for upf_file in original_upfs/*.upf; do
    # Execute the python script
    echo "Regenerating $upf_file"
    atompaw-regenerate $upf_file ./ clean
done