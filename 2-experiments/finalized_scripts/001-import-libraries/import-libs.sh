#!/bin/env bash

# To run this script, need to pip install `sssp-verify-scripts` first.
# Because it is where `import_upf_group` command is provided.
# See README.md on how to add a new pseudo library as a UPF group for verification.

for lib in \
nc-dojo-v0.4.1-std \
nc-dojo-v0.4.1-str \
nc-dojo-v0.5.0-std \
nc-dojo-v0.5.0-std-trim \
nc-psl-v1.0.0 \
nc-sg15-oncvpsp4 \
nc-spms-oncvpsp4 \
paw-actinides-marburg \
paw-jth-v1.1-std \
paw-jth-v1.1-str \
paw-lanthanides-wentzcovitch \
paw-psl-v0.x \
paw-psl-v1.0.0-high \
paw-psl-v1.0.0-low \
us-gbrv-v1.x-upf1-dont-use-this-lib-use-upf2 \
us-gbrv-v1.x-upf2 \
us-psl-v0.x \
us-psl-v1.0.0-high \
us-psl-v1.0.0-low
do
    import_upf_group --base-folder ~/project/sssp-project/sssp-verify-scripts --lib-name "$lib"
done
