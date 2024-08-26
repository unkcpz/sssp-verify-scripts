#!/bin/env python

import json
from sssp_verify_scripts.controllers import TransferabilityEOSGroupSubmissionController
from pydantic import ValidationError
from rich import print
import argparse

def launch(protocol, curate_type, library, concurrent, unit_num_cpus, ecutwfc=None, cutoff_mapping=None):
    property = "eos"
    computer = 'eiger-hq'

    # XXX: should accept group that record the mapping from pp name -> cutoff results
    # Use as cutoff_mapping and pass to controller for dynamical ACWF verification.

    target_upf_lib = f"validate/upf/candidate/{library}"

    lib_name = target_upf_lib.split('/')[-1]

    comment = f"transferability@{ecutwfc} on {property} for {curate_type} library"

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

    print("\nLaunching Convergence EOS controller ---\n")
    inputs = {
        "parent_group_label": target_upf_lib,
        "max_concurrent": concurrent,
        "pw_code": f"pw-7.2@{computer}",
        "curate_type": curate_type,
        "protocol": protocol,
        "unit_num_cpus": unit_num_cpus,
        "unit_memory_mb": unit_memory_mb,
        "unit_npool": unit_npool,
        "clean_workdir": True,
    }

    if ecutwfc is not None:
        inputs['ecutwfc'] = ecutwfc
        target_group_label = f"validate/{lib_name}/transferability/{property}/{protocol}-{curate_type}-{ecutwfc}"
        inputs["group_label"]= target_group_label
    
    if cutoff_mapping is not None:
        # XXX:flasky
        if 'efficiency' in cutoff_mapping.name:
            target_group_label = f"validate/{lib_name}/transferability/{property}/{protocol}-{curate_type}-efficiency"
        elif 'precision' in cutoff_mapping.name:
            target_group_label = f"validate/{lib_name}/transferability/{property}/{protocol}-{curate_type}-precision"
        else:
            raise ValueError(f"unexpected file: {cutoff_mapping.name}")

        cutoff_mapping = json.load(cutoff_mapping)
        inputs['cutoff_mapping'] = cutoff_mapping
        inputs["group_label"]= target_group_label

    try:
        controller = TransferabilityEOSGroupSubmissionController(
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
    parser.add_argument('--ecutwfc', type=int, default=200)
    parser.add_argument('--cutoff-mapping', type=argparse.FileType('r', encoding='utf-8'))
    parser.add_argument('--protocol', default='standard')
    parser.add_argument('--curate-type', help="`nc` or `sssp` for oxygen PP")
    # parser.add_argument('--computer', help="computer label name")
    parser.add_argument('--library', help="name of library")
    parser.add_argument('--concurrent', type=int, default=1, help="max number of concurrent")
    parser.add_argument('--ncpus', type=int, default=16, help="Number of CPUs per calculation")
    parser.add_argument('--dry-run', action='store_true', default=False)

    args = parser.parse_args()

    inputs = {
        'protocol': args.protocol,
        'curate_type': args.curate_type,
        'unit_num_cpus': args.ncpus,
        'library': args.library,
        'concurrent': args.concurrent,
    }
    
    if args.cutoff_mapping is not None:
        inputs['cutoff_mapping'] = args.cutoff_mapping
    else:
        inputs['ecutwfc'] = args.ecutwfc

    if args.dry_run:
        print(inputs)
    else:
        launch(
            **inputs,
        )


