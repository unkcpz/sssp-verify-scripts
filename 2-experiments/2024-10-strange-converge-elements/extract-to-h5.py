#!/bin/env python

from aiida import orm
import sys
import aiida
from aiida_sssp_workflow.workflows.convergence.report import ConvergenceReport
import h5py
import numpy as np
import argparse
from tqdm import tqdm
import typing as t

from aiida_sssp_workflow.utils.pseudo import extract_pseudo_info_from_filename
from aiida_sssp_workflow.workflows.convergence.bands import compute_xy as compute_xy_bands_distance 
from aiida_sssp_workflow.workflows.convergence.cohesive_energy import compute_xy as compute_xy_cohesive_energy
from aiida_sssp_workflow.workflows.convergence.eos import compute_xy as compute_xy_eos
from aiida_sssp_workflow.workflows.convergence.pressure import compute_xy as compute_xy_pressure
from aiida_sssp_workflow.workflows.transferability.eos import extract_eos as transferability_extract_eos

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def extract_convergence(property, f_compute_xy, override=False):

    with h5py.File('pp_verify_convergence.h5', 'a') as f:
        group_name = f"validate/convergence-recheck/convergence/{property}/standard/gs"
        group = orm.load_group(group_name)

        print(f"--- In processing group {group_name} ---")

        for node in tqdm(group.nodes, file=sys.stdout):
            # key is the filename of the pseudo
            filename = node.inputs.pseudo.filename
            md5 = node.inputs.pseudo.md5

            if node.exit_status != 0:
                eprint(f"node {node.uuid} verify on {filename} on/ property {property} not okay")
                continue


            # print(f"extracting {filename}")
            
            if filename in f:
                pp_grp = f[filename]
            else:
                pp_grp = f.create_group(filename)

            # Add attr:
            # - md5
            # - element
            # - functional
            # - z_valence
            # - pp_type
            # - TODO: pp_code
            # - TODO: lib_version
            # - TODO: lib_label (dojo, psl, jth etc) 
            pseudo_info = extract_pseudo_info_from_filename(filename)
            print(filename)
            parts = filename.split(".")
            if parts[5] == "jth":
                lib_name = "paw-jth-v1.1-std"
            elif "gbrv" in parts[5]:
                lib_name = "us-gbrv-v1.x-upf2"
            else:
                lib_name = "paw-psl-v0.x"

            pp_grp.attrs['element'] = pseudo_info.element
            pp_grp.attrs['functional'] = pseudo_info.functional
            pp_grp.attrs['z_valence'] = pseudo_info.z_valence
            pp_grp.attrs['pp_type'] = pseudo_info.type
            pp_grp.attrs['lib_name'] = lib_name
            pp_grp.attrs['md5'] = md5

            subgroup_name = f'convergence_{property}'
            if subgroup_name in pp_grp: 
                if override:
                    del pp_grp[subgroup_name]
                    convergence_test_group = pp_grp.create_group(subgroup_name)
                else:
                    continue
            else:
                convergence_test_group = pp_grp.create_group(subgroup_name)

            # compute the property 
            xy_dict = f_compute_xy(node)
            for key in [k for k in xy_dict.keys() if k != 'metadata']:
                convergence_test_group.create_dataset(key, data=np.array(xy_dict[key]))

def compute_xy_phonon_frequencies(
    node: orm.Node,
) -> dict[str, t.Any]:
    """From report calculate the xy data, xs are cutoffs and ys are phonon frequencies diff from reference"""
    import numpy as np

    report_dict = node.outputs.report.get_dict()
    report = ConvergenceReport.construct(**report_dict)

    reference_node = orm.load_node(report.reference.uuid)
    output_parameters_r: orm.Dict = reference_node.outputs.ph.output_parameters
    y_ref = output_parameters_r["dynamical_matrix_1"]["frequencies"]

    xs = []
    ys_relative_diff = []
    ys_omega_max = []
    ys_relative_max_diff = []
    ys_absolute_max_diff = []
    for node_point in report.convergence_list:
        if node_point.exit_status != 0:
            # TODO: log to a warning file for where the node is not finished_okay
            continue

        x = node_point.wavefunction_cutoff

        xs.append(x)

        node = orm.load_node(node_point.uuid)
        output_parameters_p: orm.Dict = node.outputs.ph.output_parameters

        y_p = output_parameters_p["dynamical_matrix_1"]["frequencies"]

        # calculate the diff
        diffs = np.array(y_p) - np.array(y_ref)
        weights = np.array(y_ref)

        relative_diff = np.sqrt(np.mean((diffs / weights) ** 2))

        omega_max = np.amax(y_p)
        absolute_max_diff = np.amax(diffs)
        relative_max_diff = np.amax(np.abs(diffs / weights))

        # Legacy modification required when configuration is `GS` which is not run for convergence anymore
        # Keep it here just for reference.
        # if configuration == "GS":
        #     # leftover setting from SSSP v1
        #     # Otherwise the phonon frequencies calculated at BZ boundary qpoint (1/2, 1/2, 1/2) are not converged.
        #     if element == "N" or element == "Cl":
        #         start_idx = 12
        #     elif element == "H" or element == "I":
        #         start_idx = 4
        #     elif element == "O":
        #         start_idx = 6
        #     else:
        #         start_idx = 0
        # else:
        #     start_idx = 0
        #
        ys_relative_diff.append(relative_diff)
        ys_omega_max.append(omega_max)
        ys_relative_max_diff.append(relative_max_diff)
        ys_absolute_max_diff.append(absolute_max_diff)

    return {
        "xs": xs,
        "ys": ys_relative_diff,
        "ys_relative_diff": ys_relative_diff,
        "ys_omega_max": ys_omega_max,
        "ys_absolute_max_diff": ys_absolute_max_diff,
        "ys_relative_max_diff": ys_relative_max_diff,
        "metadata": {
            "unit_default": "%",
        },
    }



if __name__ == "__main__":
    aiida.load_profile()

    parser = argparse.ArgumentParser()
    parser.add_argument('--override', action='store_true', default=False)

    args = parser.parse_args()

    for property in [
        # 'bands',
        # 'eos',
        'phonon_frequencies',
        # 'pressure',
        # 'cohesive_energy',
    ]:
                 
        match property:
            case 'bands':
                _compute_xy = compute_xy_bands_distance
            case 'cohesive_energy':
                _compute_xy = compute_xy_cohesive_energy
            case 'eos':
                _compute_xy = compute_xy_eos
            case 'phonon_frequencies':
                _compute_xy = compute_xy_phonon_frequencies
            case 'pressure':
                _compute_xy = compute_xy_pressure
            case _:
                raise ValueError(f"Unknow property: {property}")
    
        extract_convergence(property, _compute_xy, override=args.override)

