#!/usr/bin/env bash

for f in *.inp; do
    # Replace dft='PBE' with dft='PBESOL'
    sed -i -e "s/dft='PBE'/dft='PBESOL'/g" -e "s/dft = 'pbe'/dft='PBESOL'/g" "$f"

    # Replace file_pseudopw='...' with the input filename
    base="${f%.inp}"
    # Here we assume the original UPF filename had .upf extension
    # You can also modify to replace .pbe → .pbesol if needed
    sed -i "s|file_pseudopw='[^']*'|file_pseudopw='${base}.upf'|g" "$f"

    echo "Updated $f: PBESOL and file_pseudopw"
done
