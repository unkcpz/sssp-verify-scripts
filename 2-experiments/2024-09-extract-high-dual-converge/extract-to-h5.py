#!/bin/env python

from aiida import orm
import sys
import aiida
import h5py
import numpy as np
import argparse
from tqdm import tqdm

from aiida_sssp_workflow.utils.pseudo import extract_pseudo_info_from_filename
from aiida_sssp_workflow.workflows.convergence.bands import compute_xy as compute_xy_bands_distance 
from aiida_sssp_workflow.workflows.convergence.cohesive_energy import compute_xy as compute_xy_cohesive_energy
from aiida_sssp_workflow.workflows.convergence.eos import compute_xy as compute_xy_eos
from aiida_sssp_workflow.workflows.convergence.phonon_frequencies import compute_xy as compute_xy_phonon_frequencies
from aiida_sssp_workflow.workflows.convergence.pressure import compute_xy as compute_xy_pressure
from aiida_sssp_workflow.workflows.transferability.eos import extract_eos as transferability_extract_eos

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def extract_convergence(property, conf_dual: str, f_compute_xy, override=False):

    conf_dual_name = conf_dual.replace('/', '_')
    with h5py.File(f'pp_verify_convergence_{conf_dual_name}.h5', 'a') as f:
        group_name = f"validate/high-dual-elements/convergence/{property}/standard/{conf_dual}"
        group = orm.load_group(group_name)

        print(f"--- In processing group {group_name} ---")

        for node in tqdm(group.nodes, file=sys.stdout):
            # key is the filename of the pseudo
            filename = node.inputs.pseudo.filename
            md5 = node.inputs.pseudo.md5

            if node.exit_status != 0:
                eprint(f"node {node.uuid} verify on {filename} on/ conf {conf_dual}/ property {property} not okay")
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


if __name__ == "__main__":
    aiida.load_profile()

    parser = argparse.ArgumentParser()
    parser.add_argument('--override', action='store_true', default=False)

    args = parser.parse_args()

    for conf_dual in [
        "bcc/dual8",
        "bcc/dual12",
        "bcc/dual18",
        "dc/dual8",
        "dc/dual12",
        "dc/dual18",
    ]:
        for property in [
            'bands',
            'eos',
            'phonon_frequencies',
            'pressure',
            'cohesive_energy',
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
        
            extract_convergence(property, conf_dual, _compute_xy, override=args.override)
