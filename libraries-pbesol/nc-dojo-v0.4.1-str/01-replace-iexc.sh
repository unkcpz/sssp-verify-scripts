#!/usr/bin/env bash
set -e

for inp_file in *.inp; do
    [[ -f "$inp_file" ]] || continue

    tmp_file="$(mktemp)"

    awk '
        # Match atomic data lines (symbol + numeric fields)
        /^[A-Za-z]+\s+[0-9]/ && NF >= 6 {
            $5 = "-116133"   # replace iexc
            printf "%-3s %6.2f %4d %4d %9d %7s\n", $1, $2, $3, $4, $5, $6
            next
        }
        { print }
    ' "$inp_file" > "$tmp_file"

    mv "$tmp_file" "$inp_file"
    echo "Updated iexc in $inp_file"
done

