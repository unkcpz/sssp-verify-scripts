#!/usr/bin/env bash

for f in *.pbesol.*.upf; do
    sed -i 's/functional="SLA+PW+PSX+PSC+PBEsol"/functional="PBESOL"/' "$f"
done
