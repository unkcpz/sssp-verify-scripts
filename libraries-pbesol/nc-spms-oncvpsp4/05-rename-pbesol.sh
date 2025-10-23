#!/usr/bin/env bash

for f in *.pbe.*.upf; do
    mv "$f" "${f/.pbe./.pbesol.}"
done
