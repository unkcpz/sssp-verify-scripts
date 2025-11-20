#!/bin/env python

# -> this script might needed by Edan

# dependency:
# - "./conf_mapping.json": to map element to its "most realistic" config among bcc/fcc/sc/dc

# cp bands workchain with the target configuration to a single group
# For each pseudo lib, there are sc/dc/bcc/fcc bands groups contain bands workchain.
# Using mapping, to get the target dc/sc/bcc/fcc and put them in a single group 

# Run this script:
# `cda sssp-project` on thanos
# `verdi run 001-extract-to-a-single-group.py`

import aiida
import json

aiida.load_profile("full-24-07")

libs = [
    "nc-dojo-v0.4.1-std",
    "nc-dojo-v0.4.1-str",
    "nc-dojo-v0.5.0-std",
    "nc-spms-oncvpsp4",
    "nc-sg15-oncvpsp4",
    "us-gbrv-v1.x-upf2",
    "us-psl-v1.0.0-high",
    "us-psl-v1.0.0-low",
    "us-psl-v0.x",
    "paw-psl-v0.x",
    "paw-psl-v1.0.0-high",
    "paw-psl-v1.0.0-low",
    "paw-jth-v1.1-std",
    "paw-jth-v1.1-str",
    "paw-lanthanides-wentzcovitch",
    "paw-actinides-marburg",
]

group_label = "all-pseudos-bands-distance"
dst_group, created = Group.collection.get_or_create(label=group_label)

if created:
    print(f"Created group: {group_label}")
else:
    print(f"Group already exists: {group_label}")

nodes_to_add = []

with open("conf_mapping.json", "r") as f: 
    conf_mapping = json.load(f)

for lib in libs:
    for conf in ["dc", "sc", "fcc", "bcc"]:
        group_name = f"validate/{lib}/convergence/bands/standard/{conf}/2"
        src_group = Group.collection.get(label=group_name)

        for node in src_group.nodes:
            pp_label = node.label
            element = pp_label.split(".")[0]
            if element in ["Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]:
                continue
            if conf_mapping[element] == conf.upper():
                print(f"{pp_label}")
                nodes_to_add.append(node)

dst_group.add_nodes(nodes_to_add)

