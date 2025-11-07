#!/usr/bin/env bash
set -e

for inp_file in *.inp; do
    [[ -f "$inp_file" ]] || continue

    tmp_file="$(mktemp)"

    awk '
        # Match lines like: Symbol Number ... at least 6 fields
        /^[[:space:]]*[A-Za-z]+[[:space:]]+[0-9.]/ && NF >= 6 {
            $5 = "-116133"
        }
        { print }
    ' "$inp_file" > "$tmp_file"

    mv "$tmp_file" "$inp_file"
    echo "Updated iexc in $inp_file"
done

