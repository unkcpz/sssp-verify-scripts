#!/usr/bin/env bash

for f in *.inp; do
    sed -i "s|XC_GGA_X_PBE+XC_GGA_C_PBE|XC_GGA_X_PBE_SOL+XC_GGA_C_PBE_SOL|" "$f"

    echo "Updated $f"
done
