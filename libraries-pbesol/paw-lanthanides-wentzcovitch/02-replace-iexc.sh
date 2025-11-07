#!/usr/bin/env bash

for f in *.inp; do
    sed -i "s|GGA-PBE|XC_GGA_X_PBE_SOL+XC_GGA_C_PBE_SOL|" "$f"

    echo "Updated $f"
done
