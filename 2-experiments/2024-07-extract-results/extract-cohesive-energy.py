#!/bin/env python

from aiida_sssp_workflow.utils.pseudo import DualType, get_dual_type
import h5py
import json
from tqdm import tqdm
import pandas as pd
import sys
from tabulate import tabulate
from aiida_sssp_workflow.utils.protocol import get_protocol

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Load the dataset of convergence results
converge_h5 = h5py.File('./pp_verify_convergence.h5')

# traverse once to collect mapping of lib -> all its PPs
lib_pps_mapping = {}

def curated_by_lib(name: str, obj):
    # only get result for first layer
    if '/' in name:
        return
    lib_name = obj.attrs.get('lib_name')
    if lib_name is None:
        raise ValueError(f"lib_name attr of {obj} is None")

    lib_pps_mapping.setdefault(lib_name, []).append(name)


converge_h5.visititems(curated_by_lib)

def extract_cohesive_energy_at_max(lib_name):
    pps = lib_pps_mapping[lib_name]

    df = pd.DataFrame(columns=properties + ['max'])

    mapping = {}
    for pp_name in tqdm(pps, file=sys.stdout):
        dataset = converge_h5[pp_name]

        md5 = dataset.attrs.get('md5')
        element = dataset.attrs.get('element')

        if md5 is None:
            raise ValueError(f"md5 of {dataset} is None")

        # print(f"------> Pseudopotential = {pp_name}")

        property = "cohesive_energy"

        try:
            xs = dataset[f'convergence_{property}']['xs'][()]
            ys = dataset[f'convergence_{property}']['ys'][()]
        except KeyError:
            eprint(f"not able to get {property} of {pp_name}")
            cutoff = None
        else:
            if xs[-1] == 200:
                y = ys[-1]
            else:
                continue
        
        yys.append(y)

        df.loc[pp_name] = cutoffs + [max_cutoff]
        mapping[pp_name] = {
            'md5': md5,
            'element': element,
            'cohesive_energy_at_max': y,
        }

    table = tabulate(df, headers='keys', tablefmt='pretty')
    with open(f'cohesive_{lib_name}-{protocol}.txt', 'w') as f:
        f.write(table)

    with open(f'cohesive_{lib_name}-{protocol}.json', 'w', encoding='utf-8') as f:
        json.dump(mapping, f, indent=4)

if __name__ == '__main__':

    for lib_name in [
        "nc-dojo-v0.4.1-std",
        "paw-jth-v1.1-std",
        "us-gbrv-v1.x-upf2",
        "us-psl-v1.0.0-high",
        "nc-spms-oncvpsp4",
    ]:
        extract(lib_name)
