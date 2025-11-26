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
converge_h5_dc = h5py.File('./pp_verify_convergence_dc.h5')
converge_h5_bcc = h5py.File('./pp_verify_convergence_bcc.h5')
converge_h5_fcc = h5py.File('./pp_verify_convergence_fcc.h5')

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

converge_h5_dc.visititems(curated_by_lib)
converge_h5_bcc.visititems(curated_by_lib)
converge_h5_fcc.visititems(curated_by_lib)

if __name__ == '__main__':

    mapping = {} # {'Au': {'fcc': xxx, 'bcc': xxx}}

    for conf, converge_h5 in [('bcc', converge_h5_bcc), ('fcc', converge_h5_fcc), ('dc', converge_h5_dc)]:
        for lib_name in [
                "nc-dojo-v0.4.1-std",
                "nc-dojo-v0.4.1-str",
                "nc-dojo-v0.5.0-std",
                "nc-sg15-oncvpsp4",
                "nc-spms-oncvpsp4",
                "paw-jth-v1.1-std",
                "paw-jth-v1.1-str",
                "paw-lanthanides-wentzcovitch",
                "paw-psl-v0.x",
                "paw-psl-v1.0.0-high",
                "paw-psl-v1.0.0-low",
                "paw-actinides-marburg",
                "us-gbrv-v1.x-upf2",
                "us-psl-v1.0.0-high",
                "us-psl-v1.0.0-low",
                "us-psl-v0.x",
        ]:
            pps = lib_pps_mapping[lib_name]

            df = pd.DataFrame(columns=['element', 'cohesive_energy'])

            for pp_name in tqdm(pps, file=sys.stdout):
                try:
                    dataset = converge_h5[pp_name]
                except:
                    continue

                md5 = dataset.attrs.get('md5')
                element = dataset.attrs.get('element')

                if md5 is None:
                    raise ValueError(f"md5 of {dataset} is None")

                # print(f"------> Pseudopotential = {pp_name}")

                property = "cohesive_energy"

                try:
                    xs = dataset[f'convergence_{property}']['xs'][()]
                    ys = dataset[f'convergence_{property}']['ys_cohesive_energy_per_atom'][()]
                except KeyError:
                    eprint(f"not able to get {property} of {pp_name}")
                    y = None
                else:
                    if xs[-1] == 200:
                        y = ys[-1]
                        print(element)
                        print(y)
                        if element not in mapping:
                            mapping[element] = {}

                        try:
                            current_stable = mapping[element][conf]
                        except KeyError:
                            mapping[element][conf] = y
                        else:
                            if y < current_stable:
                                mapping[element][conf] = y
                    else:
                        y = None

    print(mapping)
