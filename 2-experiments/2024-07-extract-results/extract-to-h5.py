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

def extract_convergence(lib_name, property, f_compute_xy, override=False):

    with h5py.File('pp_verify_convergence.h5', 'a') as f:
        group_name = f"validate/{lib_name}/convergence/{property}/standard/dc"
        group = orm.load_group(group_name)

        print(f"--- In processing group {group_name} ---")

        for node in tqdm(group.nodes, file=sys.stdout):
            # key is the filename of the pseudo
            filename = node.inputs.pseudo.filename
            md5 = node.inputs.pseudo.md5

            if node.exit_status != 0:
                eprint(f"node {node.uuid} verify on {filename} not okay")
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

def extract_eos(lib_name, override=False):

    criteria = 200 # XXX: 200, eff (recommended cutoff using eff criteria), prec (..)
    with h5py.File(f'pp_verify_transferability_eos_{criteria}.h5', 'a') as f:
        group_name = f"validate/{lib_name}/transferability/eos/standard-sssp-{criteria}"
        group = orm.load_group(group_name)

        print(f"--- In processing group {group_name} ---")

        for node in tqdm(group.nodes, file=sys.stdout):
            # key is the filename of the pseudo
            filename = node.inputs.pseudo.filename
            md5 = node.inputs.pseudo.md5

            if node.exit_status != 0:
                # XXX: this should loose a bit since any configuraiton failed will result to 
                # no data for other success confs
                eprint(f"node {node.uuid} verify on {filename} not okay")
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
            pp_grp.attrs['element'] = pseudo_info.element
            pp_grp.attrs['functional'] = pseudo_info.functional
            pp_grp.attrs['z_valence'] = pseudo_info.z_valence
            pp_grp.attrs['pp_type'] = pseudo_info.type
            pp_grp.attrs['lib_name'] = lib_name
            pp_grp.attrs['md5'] = md5

            subgroup_name = f'transferability_eos'
            if subgroup_name in pp_grp: 
                if override:
                    del pp_grp[subgroup_name]
                    transferability_eos_group = pp_grp.create_group(subgroup_name)
                else:
                    continue
            else:
                transferability_eos_group = pp_grp.create_group(subgroup_name)

            # Store dataset in every group of a configuration
            eos_dict, metric_dict = transferability_extract_eos(node)
            for conf in [k for k in eos_dict.keys() if k != 'metadata']:
                c_group = transferability_eos_group.create_group(conf)

                volume_energy_dict = eos_dict[conf]
                volumes = volume_energy_dict['volumes']
                energies = volume_energy_dict['energies']
                c_group.create_dataset('volumes', data=volumes)
                c_group.create_dataset('energies', data=energies)

                c_group.attrs['delta'] = metric_dict[conf]['delta/natoms']
                c_group.attrs['delta1'] = metric_dict[conf]['delta1']
                c_group.attrs['nu'] = metric_dict[conf]['rel_errors_vec_length']


if __name__ == "__main__":
    aiida.load_profile()

    parser = argparse.ArgumentParser()
    parser.add_argument('--override', action='store_true', default=False)
    parser.add_argument('--group', help='`convergence` or `eos`')

    args = parser.parse_args()

    for lib_name in [
        "nc-dojo-v0.4.1-std",
        "paw-jth-v1.1-std",
        "us-gbrv-v1.x-upf2",
        "us-psl-v1.0.0-high",
        "nc-spms-oncvpsp4",
    ]:
        if args.group == 'convergence':
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

                extract_convergence(lib_name, property, _compute_xy, override=args.override)
        elif args.group == 'eos':
            extract_eos(lib_name, override=args.override)
        else:
            raise ValueError(f"unknow group flag {args.group}")
