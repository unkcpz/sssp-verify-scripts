#!/bin/env python

from copy import deepcopy
import h5py
import json
import sys
from aiida_sssp_workflow.utils.element import ALL_ELEMENTS
from aiida_sssp_workflow.calculations.calculate_metric import rel_errors_vec_length

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


lib_abbr_name_mapping = {
    'nc-dojo-v0.4.1-std': 'DOJO-041-std',
    'nc-spms-oncvpsp4': 'SPMS',
    "nc-dojo-v0.4.1-str": 'DOJO-041-str',
    "nc-dojo-v0.5.0-std": 'DOJO-050-std',
    "nc-sg15-oncvpsp4": 'SG15',
    'us-gbrv-v1.x-upf2': 'GBRV-1.X',
    'us-psl-v1.0.0-high': 'PSL-US-v1-high',
    "us-psl-v1.0.0-low": 'PSL-US-v1-low',
    "us-psl-v0.x": 'PSL-US-v0.x',
    'paw-jth-v1.1-std': 'JTH-1.1-std',
    "paw-jth-v1.1-str": 'JTH-1.1-str',
    "paw-lanthanides-wentzcovitch": 'Wentzcovitch',
    "paw-psl-v0.x": 'PSL-PAW-v0.x',
    "paw-psl-v1.0.0-high": 'PSL-PAW-v1-high',
    "paw-psl-v1.0.0-low": 'PSL-PAW-v1-low',
    "paw-actinides-marburg": 'MARBURG',
}

REF_DATA = {}
with open("./results-oxides-verification-PBE-v1-AE-average.json", "r") as fh:
    data = json.load(fh)
    data = data["BM_fit_data"]

    for element in ALL_ELEMENTS:
        if element not in REF_DATA:
            REF_DATA[element] = {}

        for c in ['XO', 'XO2', 'XO3', 'X2O', 'X2O3', 'X2O5']:
            try:
                d = data[f"{element}-{c}"]
            except KeyError:
                continue
            REF_DATA[element][c] = {
                'V0': d['min_volume'],
                'B0': d['bulk_modulus_ev_ang3'],
                'B1': d['bulk_deriv'],
            }
            
with open("./results-unaries-verification-PBE-v1-AE-average.json", "r") as fh:
    data = json.load(fh)
    data = data["BM_fit_data"]

    for element in ALL_ELEMENTS:
        if element not in REF_DATA:
            REF_DATA[element] = {}

        for c in ['BCC', 'FCC', 'SC', 'DC']:
            if c == 'DC':
                cc = 'Diamond'
            else:
                cc = c

            try:
                d = data[f"{element}-X/{cc}"]
            except KeyError:
                continue
            REF_DATA[element][c] = {
                'V0': d['min_volume'],
                'B0': d['bulk_modulus_ev_ang3'],
                'B1': d['bulk_deriv'],
            }

def extract(element, element_pps_mapping) -> dict:

    pps_info = {}

    try:
        pps = element_pps_mapping[element]
    except KeyError:
        return {}

    pps_info['REF'] = REF_DATA[element]

    if element == "O":
        pps_info['REF']["XO"]["V0"] *= 2
        pps_info['REF']["X2O5"]["V0"] *= 2

    for pp_name in pps:
        print(f"------> Pseudopotential = {pp_name}")

        eos_dataset = eos_h5[pp_name]
        lib_name = eos_dataset.attrs.get('lib_name')
        lib_name_abbr = lib_abbr_name_mapping[lib_name]

        pp_info = {}
        for conf, data in eos_dataset['transferability_eos'].items():
            volumes = data['volumes']
            energies = data['energies']

            V0 = data.attrs.get('V0') 
            B0 = data.attrs.get('B0') 
            B1 = data.attrs.get('B1') 
            E0 = data.attrs.get('E0') 

            # acwf ref use premitive cell, I use same cell as other elements.
            if element == "O" and conf in ["XO", "X2O5"]:
                # acwf ref use primitive cell for XO and X2O5 the cell size is half
                ref_V0 = pps_info['REF'][conf]["V0"]
                ref_B0 = pps_info['REF'][conf]["B0"]
                ref_B1 = pps_info['REF'][conf]["B1"]
                nu = rel_errors_vec_length(ref_V0, ref_B0, ref_B1, V0, B0, B1)
            else:
                nu = data.attrs.get('nu')

            pp_info[conf] = {
                'volumes': list(volumes),
                'energies': list(energies),
                'nu': nu,
                'V0': V0, 
                'B0': B0, 
                'B1': B1, 
                'E0': E0, 
            }

        pps_info[lib_name_abbr] = pp_info

    return pps_info
            
if __name__ == "__main__":
    # traverse once to collect mapping of element -> all PPs
    eos_h5 = h5py.File('./pp_verify_transferability_eos_200.h5')
    element_pps_mapping = {}

    def curated_by_element(name: str, obj):
        # only get result for first layer
        if '/' in name:
            return
        element = obj.attrs.get('element')
        if element is None:
            raise ValueError(f"element attr of {obj} is None")

        element_pps_mapping.setdefault(element, []).append(name)

    eos_h5.visititems(curated_by_element)

    with open("eos.json", "w") as fh:
        info = {}
        for element in ALL_ELEMENTS:
        # for element in ["Fe"]:
            epp_info = extract(element, element_pps_mapping)
            
            info[element] = epp_info

        json.dump(info, fh, indent=4)
