#!/bin/env python

from aiida_sssp_workflow.utils.pseudo import DualType, get_dual_type
import h5py
import json
from tqdm import tqdm
import pandas as pd
import sys
from tabulate import tabulate
from aiida_sssp_workflow.utils.protocol import get_protocol

def compute_recommended_cutoffs(xs, ys, criteria) -> int:
    for x, y in zip(reversed(xs), reversed(ys)):
        cutoff = x
        if y > criteria:
            break

    return cutoff

def get_criteria(protocol, property):
    criteria = get_protocol(category='criteria', name=protocol)

    return criteria[property]['bounds'][1]


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

def extract(lib_name, protocol):
    pps = lib_pps_mapping[lib_name]

    properties = ['bands', 'eos', 'phonon_frequencies', 'pressure', 'cohesive_energy']
    df = pd.DataFrame(columns=properties + ['max'])

    cutoff_mapping = {}
    for pp_name in tqdm(pps, file=sys.stdout):
        dataset = converge_h5[pp_name]

        md5 = dataset.attrs.get('md5')
        element = dataset.attrs.get('element')
        pp_type = dataset.attrs.get('pp_type')
        match get_dual_type(pp_type, element):
            case DualType.NC:
                dual = 4
            case DualType.AUGLOW:
                dual = 8
            case DualType.AUGHIGH:
                dual = 18

        if md5 is None:
            raise ValueError(f"md5 of {dataset} is None")

        # print(f"------> Pseudopotential = {pp_name}")

        cutoffs = []
        for property in properties:
            criteria = get_criteria(protocol, property)

            try:
                xs = dataset[f'convergence_{property}']['xs'][()]
                ys = dataset[f'convergence_{property}']['ys'][()]
            except KeyError:
                eprint(f"not able to get {property} of {pp_name}")
                cutoff = None
            else:
                cutoff = compute_recommended_cutoffs(xs, ys, criteria)
            
            cutoffs.append(cutoff)

        if None in cutoffs:
            max_cutoff = 200
        else:
            max_cutoff = max(cutoffs)

        df.loc[pp_name] = cutoffs + [max_cutoff]
        cutoff_mapping[pp_name] = {
            'md5': md5,
            'cutoffs': (int(max_cutoff), int(max_cutoff * dual)),
        }

    table = tabulate(df, headers='keys', tablefmt='pretty')
    with open(f'{lib_name}-{protocol}.txt', 'w') as f:
        f.write(table)

    with open(f'{lib_name}-{protocol}.json', 'w', encoding='utf-8') as f:
        json.dump(cutoff_mapping, f, indent=4)

if __name__ == '__main__':

    for protocol in ['efficiency', 'precision']:
        for lib_name in [
            "nc-dojo-v0.4.1-std",
            "paw-jth-v1.1-std",
            "us-gbrv-v1.x-upf2",
            "us-psl-v1.0.0-high",
            "nc-spms-oncvpsp4",
        ]:
            extract(lib_name, protocol)
