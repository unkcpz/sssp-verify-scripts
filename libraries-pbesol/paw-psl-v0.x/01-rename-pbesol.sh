#!/usr/bin/env bash

for f in *.pbe.*.inp; do
    mv "$f" "${f/.pbe./.pbesol.}"
done
