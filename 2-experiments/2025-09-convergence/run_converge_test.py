#!/bin/env python

import argparse
from aiida_sssp_workflow.utils import get_protocol
from pydantic import ValidationError

from sssp_verify_scripts.controllers.rho_cohesive_energy import (
    ConvergenceCohesiveEnergyGroupSubmissionController,
)
from sssp_verify_scripts.controllers.rho_bands import (
    ConvergenceBandsGroupSubmissionController,
)
from sssp_verify_scripts.controllers.rho_eos import (
    ConvergenceEOSGroupSubmissionController,
)
from sssp_verify_scripts.controllers.rho_phonon_frequencies import (
    ConvergencePhononFrequenciesGroupSubmissionController,
)
from sssp_verify_scripts.controllers.rho_pressure import (
    ConvergencePressureGroupSubmissionController,
)


# batch=2: start at 2024-12-03
def launch(
    protocol,
    property,
    configuration,
    experiment,
    library,
    wavefunction_cutoff_mapping,
    concurrent,
    unit_num_cpus,
    batch=2,
):
    target_upf_lib = f"validate/upf/candidate/{library}"

    lib_name = target_upf_lib.split("/")[-1]
    target_group_label = f"validate/{lib_name}/convergence/{property}/{protocol}/{configuration.lower()}/{batch}"

    comment = f"convergence test@ {property}, with {configuration}. experiment mode? ({experiment})"

    print(f"{comment}")

    if unit_num_cpus % 4 != 0:
        raise ValueError(f"Expect times of 4: got {unit_num_cpus}")

    mem_per_cpu = 10000  # mb
    unit_memory_mb = mem_per_cpu * unit_num_cpus  # mb
    unit_npool = unit_num_cpus / 4

    control_mapping = {
        "standard": "standard",
        "experiment": "experiment",
    }

    # XXX: convergence protocol should have a independent CLI option to override
    convergence_mapping = {
        "standard": "balanced",
        "experiment": "balanced",
    }

    controller_inputs = {
        "group_label": target_group_label,
        "parent_group_label": target_upf_lib,
        "max_concurrent": concurrent,
        "pw_code": "qe-7.3-pw-gf@thor",
        "protocol": convergence_mapping[protocol],
        "configuration": configuration,
        "element_wavefunction_cutoff_mapping": wavefunction_cutoff_mapping,
        "unit_num_cpus": unit_num_cpus,
        "unit_memory_mb": unit_memory_mb,
        "unit_npool": unit_npool,
        "clean_workdir": True,
    }

    if property == "cohesive_energy":
        _SubmissionController = ConvergenceCohesiveEnergyGroupSubmissionController
    elif property == "bands":
        _SubmissionController = ConvergenceBandsGroupSubmissionController
    elif property == "phonon_frequencies":
        _SubmissionController = ConvergencePhononFrequenciesGroupSubmissionController
        controller_inputs["ph_code"] = "qe-7.3-ph-gf@thor"
    elif property == "pressure":
        _SubmissionController = ConvergencePressureGroupSubmissionController
    elif property == "eos":
        _SubmissionController = ConvergenceEOSGroupSubmissionController
    else:
        raise ValueError(f"Unknow property {property}")

    print(f"\nLaunching Convergence {property.upper()} controller ---\n")
    try:
        controller = _SubmissionController(
            **controller_inputs,
        )
    except ValidationError as exc:
        print(repr(exc.errors()))
    else:
        controller.submit_new_batch(verbose=True)


if __name__ == "__main__":
    import aiida
    import json

    aiida.load_profile()

    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", default="standard")
    parser.add_argument("--property")
    parser.add_argument(
        "--experiment",
        action="store_true",
        default=False,
        help="Using dense cutoff grid. default=False",
    )
    parser.add_argument(
        "--configuration",
        default="DC",
        help="configuration to run convergence, default=DC",
    )
    # parser.add_argument('--computer', help="computer label name")
    parser.add_argument("--library", help="name of library")
    parser.add_argument("--library-ecut-mapping", help="ecut mapping file")
    parser.add_argument(
        "--concurrent", type=int, default=1, help="max number of concurrent"
    )
    parser.add_argument(
        "--ncpus", type=int, default=16, help="Number of CPUs per calculation"
    )
    parser.add_argument("--dry-run", action="store_true", default=False)

    args = parser.parse_args()

    with open(args.library_ecut_mapping, mode="r") as f:
        wavefunction_cutoff_mapping = json.load(f)

    if args.dry_run:
        print(args)
    else:
        launch(
            protocol=args.protocol,
            property=args.property,
            # computer=args.computer,
            library=args.library,
            wavefunction_cutoff_mapping=wavefunction_cutoff_mapping,
            unit_num_cpus=args.ncpus,
            configuration=args.configuration,
            experiment=args.experiment,
            concurrent=args.concurrent,
        )
