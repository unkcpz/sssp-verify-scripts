#!/bin/env python

# -> this script may need by Edan

# The result of this script is to generate a json file "pseudo-to-bands-uuid-mapping.json" that maps from pseduo file name into
# the uuid of the WorkChain node and Bands node.

# The json will looks like
# {
#     "Th.paw.pbe.z_12.ld1.uni-marburg.v0.upf": {
#         "banddata-node": "c16ce261-cb7c-4f63-b6d6-870e004e9910",
#         "bandpara-node": "c54d0e20-cdb7-4800-a84e-68adb40f2419",
#         "workchain-node": "0d4cb51c-e8da-4a96-8000-b0c7875a747b"
#     },
#     "Pa.paw.pbe.z_13.ld1.uni-marburg.v0.upf": {
#         "banddata-node": "d56ac149-7e74-4d2d-b415-23bd1896f51c",
#         "bandpara-node": "da269eb8-7979-468b-b9c8-329aa19f6f47",
#         "workchain-node": "03e7a8c2-c44b-43c1-b728-dc4d58e6ad4b"
#     },
#     "Np.paw.pbe.z_15.ld1.uni-marburg.v0.upf": {
#         "banddata-node": "dbb77888-513c-4b27-b0d7-3a4808381e32",
#         "bandpara-node": "eb4eb050-5800-4868-a161-d86b1889d0f6",
#         "workchain-node": "b0f2d9e5-edfa-4ab2-bad5-7295bb474cb9"
#     },
#     "U.paw.pbe.z_14.ld1.uni-marburg.v0.upf": {
#         "banddata-node": "ff278eb1-acd7-4b22-83da-65e12219925d",
#         "bandpara-node": "339877ee-670d-4526-b3c1-8d10ad55d48f",
#         "workchain-node": "40daa5e6-0e3f-4565-bed2-47ed7db1ff82"
#     }
# }

# Run this script:
# `cda sssp-project` on thanos
# `cd` to this folder
# `verdi run 002-generate-mapping-pseudo-to-uuid.py`

import aiida
import json

aiida.load_profile("full-24-07")

group_name = "all-pseudos-bands-distance"
group = Group.collection.get(label=group_name)

mapping = {}

for node in group.nodes:
    pp_label = node.label
    print(f"on {pp_label}")

    matches = [
        link.node
        for link in node.base.links.get_outgoing(node_class=WorkChainNode)
        if link.link_label.startswith("cutoffs_200_")
    ]

    band_wc = matches[0]
    band_parameters_node = band_wc.outputs["bands"]["band_parameters"]
    band_structure_node = band_wc.outputs["bands"]["band_structure"]
    mapping[pp_label] = {
        "banddata-node": str(band_structure_node.uuid),
        "bandpara-node": str(band_parameters_node.uuid),
        "workchain-node": str(band_wc.uuid),
    }

with open("pseudo-to-bands-uuid-mapping.json", mode="w") as fh:
    json.dump(mapping, fh, indent=4)
