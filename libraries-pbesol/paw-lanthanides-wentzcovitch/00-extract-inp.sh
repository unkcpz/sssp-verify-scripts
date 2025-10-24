#!/usr/bin/env bash

for f in *.upf; do
    out="${f%.upf}.inp"
    
    # Extract content after "UPF file from ATOMPAW code with following input" up to </PP_INFO>
    awk '
        /UPF file from ATOMPAW code with following input/ {flag=1; next} 
        /<\/PP_INFO>/ {flag=0; next} 
        flag {print}' "$f"  | sed 's/^[[:space:]]*//' > "$out"
    
    echo "Created $out"
done
