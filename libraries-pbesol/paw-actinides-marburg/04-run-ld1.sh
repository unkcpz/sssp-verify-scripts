#!/usr/bin/env bash
set -euo pipefail

# Create tmp folder in current directory
mkdir -p tmp
mkdir -p inps

# Loop over all .inp files in current directory
for inp in *.inp; do
    [[ -f "$inp" ]] || continue

    base="${inp%.inp}"   # e.g., Ag.inp → Ag
    echo "Processing $inp ..."

    # Move the file into tmp/
    cp "$inp" inps/
    mv "$inp" tmp/

    # Enter tmp
    cd tmp

    # Run Docker command
    docker run -i \
        -v "$PWD":/workdir:rw \
        -w /workdir \
        -u "$(id -u):$(id -g)" \
        pspgen/ld1:main \
        sh -c "ld1.x < \"$inp\" > \"${base}.out\""

    # Move the .upf file back to parent folder
    mv "${base}.upf" ..

    # Optional: move the full .out back as well
    mv "${base}.out" ..

    # Optional: remove the inp from tmp/
    rm -f "$inp"

    # Go back to parent directory
    cd ..

done

# Remove tmp folder if empty
rmdir tmp 2>/dev/null || true

echo "✅ All runs completed. Output files are in $(pwd)/"

