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

    # Run Docker command using atompaw
    if ! docker run -i \
        -v "$PWD":/workdir:rw \
        -w /workdir \
        -u "$(id -u ${USER}):$(id -g ${USER})" \
        ghcr.io/pspgen/atompaw:latest \
        sh -c "atompaw < \"$inp\" > \"${base}.out\""; then
        echo "⚠️ Docker command failed for $inp, skipping to next file."
        cd ..
        continue
    fi

    # Extract element symbol from base (everything before first dot)
    element="${base%%.*}"

    # Move the generated UPF to parent folder, renaming it
    mv "${element}.GGA-PBE-paw.UPF" "../${base}.upf"

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

