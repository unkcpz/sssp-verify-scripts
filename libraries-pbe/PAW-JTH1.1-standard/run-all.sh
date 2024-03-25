#!/bin/bash

input_folder="."

for upf_file in *.upf; do
    # Execute the python script
    echo "Regenerating $upf_file"
    python regenerate.py $upf_file inputs
done