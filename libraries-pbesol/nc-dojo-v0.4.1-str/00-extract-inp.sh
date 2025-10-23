#!/usr/bin/env bash
for f in *.upf; do
    out="${f%.upf}.inp"
    awk '/<PP_INPUTFILE>/{flag=1;next}/<\/PP_INPUTFILE>/{flag=0}flag' "$f" > "$out"
    echo "Created $out"
done
