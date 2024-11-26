#!/bin/env python

import json
from sssp_verify_scripts.controllers.magnetization import MagnetizationGroupSubmissionController
from pydantic import ValidationError
from rich import print
import argparse
import numpy as np

def launch(protocol, library, concurrent, unit_num_cpus, ecutwfc, configuration, step=1, cutoff_mapping=None):
    property = "magnetization"
    computer = 'eiger-hq'

    # XXX: should accept group that record the mapping from pp name -> cutoff results
    # Use as cutoff_mapping and pass to controller for dynamical ACWF verification.

    target_upf_lib = f"validate/upf/candidate/{library}"

    lib_name = target_upf_lib.split('/')[-1]

    comment = f"transferability@(ecutwfc={ecutwfc}) on {property} on {configuration}"

    print(f"{comment}")

    # unit_num_cpus = 128
    # unit_memory_mb = 120000 # not used
    # unit_npool = 16

    # hq
    if unit_num_cpus % 4 != 0:
        raise ValueError(f"Expect times of 4: got {unit_num_cpus}")

    mem_per_cpu = 3500 # mb
    unit_memory_mb = mem_per_cpu * unit_num_cpus # mb
    unit_npool = unit_num_cpus / 4

    # computer = 'daint-hq'
    # unit_num_cpus = 18
    # unit_memory_mb = 60000 # mb
    # unit_npool = 2
    scale_list = ((1.0 + 0.01 * np.arange(-6, 16, step)) ** 3).tolist()

    print("\nLaunching transferability magnetization controller ---\n")
    inputs = {
        "parent_group_label": target_upf_lib,
        "max_concurrent": concurrent,
        "pw_code": f"pw-7.2@{computer}",
        "scale_list": scale_list,
        "ecutwfc": ecutwfc,
        "protocol": protocol,
        "unit_num_cpus": unit_num_cpus,
        "unit_memory_mb": unit_memory_mb,
        "unit_npool": unit_npool,
        "configuration": configuration,
        "clean_workdir": True,
    }

    target_group_label = f"validate/{lib_name}/transferability/{property}/{configuration}/{protocol}-{ecutwfc}"
    inputs["group_label"]= target_group_label

    try:
        controller = MagnetizationGroupSubmissionController(
            **inputs,
        )
    except ValidationError as exc:
        print(repr(exc.errors()))
    else:
        controller.submit_new_batch(verbose=True)

if __name__ == "__main__":
    import aiida
    
    aiida.load_profile()

    parser = argparse.ArgumentParser()
    parser.add_argument('--protocol', default='standard')
    parser.add_argument('--configuration', default='FCC')
    # parser.add_argument('--computer', help="computer label name")
    parser.add_argument('--ecutwfc', type=int, default=80, help="wavefunction cutoff")
    parser.add_argument('--step', type=float, default=0.1, help="step of points in volume variant")
    parser.add_argument('--library', help="name of library")
    parser.add_argument('--concurrent', type=int, default=1, help="max number of concurrent")
    parser.add_argument('--ncpus', type=int, default=16, help="Number of CPUs per calculation")
    parser.add_argument('--dry-run', action='store_true', default=False)

    args = parser.parse_args()

    inputs = {
        'protocol': args.protocol,
        'unit_num_cpus': args.ncpus,
        'library': args.library,
        'concurrent': args.concurrent,
        'configuration': args.configuration,
        'step': args.step,
        'ecutwfc': args.ecutwfc,
    }
    
    if args.dry_run:
        print(inputs)
    else:
        launch(
            **inputs,
        )


