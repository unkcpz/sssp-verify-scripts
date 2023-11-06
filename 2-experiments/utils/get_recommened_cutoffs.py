#!/usr/bin/env runaiida
"""The script to generate the recommended cutoffs for the label of pseudos
input is a group of convergence workchains, the output is a dict mapping
from the pseudo label to the recommended cutoffs, which can be used to run the batch of 
precision measure workchains.
"""
import argparse

from pprint import pprint
from aiida import orm

# group name from arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("group", type=str, help="group name")

args = parser.parse_args()
group_name = args.group

# load group
group = orm.Group.get(label=group_name)

out = {}

# get all workchains
for wc in group.nodes:
    # get the label of the pseudo
    pseudo_label = wc.base.extras.all['label'].split(' ')[-1]
    # get the recommended cutoffs by loop over all properties
    max_ecutwfc = 0
    max_ecutrho = 0
    for conv_out in wc.outputs.convergence.values():
        output_parameters = conv_out.output_parameters.get_dict()

        ecutwfc = output_parameters['wavefunction_cutoff']
        ecturho = output_parameters['chargedensity_cutoff']

        max_ecutwfc = max(max_ecutwfc, ecutwfc)
        max_ecutrho = max(max_ecutrho, ecturho)

    recommended_cutoffs = {
        "wavefunction_cutoff": max_ecutwfc,
        "chargedensity_cutoff": max_ecutrho,
    }

    out[pseudo_label] = tuple(recommended_cutoffs.values())

pprint(out)