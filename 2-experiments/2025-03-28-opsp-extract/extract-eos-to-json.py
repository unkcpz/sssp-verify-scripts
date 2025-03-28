#!/bin/env python

from aiida import orm
import sys
import numpy as np
import argparse
import json
from aiida_sssp_workflow.utils.element import ALL_ELEMENTS
from aiida_sssp_workflow.calculations.calculate_metric import rel_errors_vec_length

from aiida_sssp_workflow.workflows.transferability.eos import (
    extract_eos as transferability_extract_eos,
)

def extract_custom_name(filename: str) -> str:
    parts = filename.split(".")

    return ".".join(parts[1:-2])

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

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

def extract_eos(element, pks):
    pps_info = {}
    pps_info['REF'] = REF_DATA[element]

    if element == "O":
        pps_info['REF']["XO"]["V0"] *= 2
        pps_info['REF']["X2O5"]["V0"] *= 2

    with open(f"{element}.EOS.json", "w") as fh:

        for pk in pks:
            node = orm.load_node(pk)
            filename = node.inputs.pseudo.filename
            if node.exit_status != 0:
                eprint(f"node {node.uuid} verify on {filename} not okay")
                raise
            
            print(f"------> Pseudopotential = {filename}")
            custom_name = extract_custom_name(filename)

            eos_dict, birch_murnaghan_fit, metric_dict = transferability_extract_eos(
                node
            )

            pp_info = {}
            for conf in [k for k in eos_dict.keys() if k != "metadata"]:
                volume_energy_dict = eos_dict[conf]
                volumes = volume_energy_dict["volumes"]
                energies = volume_energy_dict["energies"]

                read_nu = metric_dict[conf]["rel_errors_vec_length"]

                V0, B0, B1 = list(metric_dict[conf]["birch_murnaghan_results"])
                E0 = birch_murnaghan_fit[conf]["energy0"]

                # acwf ref use premitive cell, I use same cell as other elements.
                if element == "O" and conf in ["XO", "X2O5"]:
                    # acwf ref use primitive cell for XO and X2O5 the cell size is half
                    ref_V0 = REF_DATA[element][conf]["V0"]
                    ref_B0 = REF_DATA[element][conf]["B0"]
                    ref_B1 = REF_DATA[element][conf]["B1"]
                    nu = rel_errors_vec_length(ref_V0, ref_B0, ref_B1, V0, B0, B1)
                else:
                    nu = read_nu

                pp_info[conf] = {
                    'volumes': list(volumes),
                    'energies': list(energies),
                    'nu': nu,
                    'V0': V0, 
                    'B0': B0, 
                    'B1': B1, 
                    'E0': E0, 
                }

            pps_info[custom_name] = pp_info

        json.dump(pps_info, fh, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--element', default=False)
    parser.add_argument('pks', nargs='+', default=False)

    args = parser.parse_args()

    extract_eos(element=args.element, pks=args.pks)
