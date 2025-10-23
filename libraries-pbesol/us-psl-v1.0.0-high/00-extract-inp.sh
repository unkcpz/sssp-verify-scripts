#!/usr/bin/env bash

for f in *.upf; do
    out="${f%.upf}.inp"
    
    # Extract lines starting from &input up to but not including the closing ]]></PP_INPUTFILE>
    awk '/&input/{flag=1} /\]\]><\/PP_INPUTFILE>/{flag=0} flag {print}' "$f" > "$out"
    
    echo "Created $out"
done

