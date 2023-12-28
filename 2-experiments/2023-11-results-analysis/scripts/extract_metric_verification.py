"""Extract JSON from the AiiDA group of metric measure verification calculations.
"""
import tqdm
import os
import json
import numpy as np

from collections import Counter
from eos_utils.eosfit_31_adapted import BM, echarge

from aiida import load_profile, orm
from argparse import ArgumentParser

__version__ = '0.1.0'

def main():
    parser = ArgumentParser()
    parser.add_argument("libname", help="Name of the pseudo library of metric measure verification calculations.")
    parser.add_argument("setname", help="Set name of configuraitons, unaries or oxides")
    parser.add_argument("output", help="Name of the output JSON file.")
    args = parser.parse_args()

    WORKFLOWS_GROUP_LABEL = f'{args.libname}'

    # Get all nodes in the output group (EOS workflows)
    group_node_query = orm.QueryBuilder().append(
        orm.Group, filters={'label': WORKFLOWS_GROUP_LABEL}, tag='groups',
    ).append(orm.Node, project='*', with_group='groups')
    group_node_query.distinct()
    wf_nodes = group_node_query.all(flat=True)
    
    states = []
    data_to_print = {}
    warning_lines = []

    uuid_mapping = {}
    all_missing_outputs = {}
    completely_off = []
    failed_wfs = []
    all_eos_data = {}
    all_stress_data = {}
    all_BM_fit_data = {}
    num_atoms_in_sim_cell = {}

    if args.setname == "unaries":
        conf_list = ["BCC", "FCC", "SC", "Diamond"]
    elif args.setname == "oxides":
        conf_list = ["XO", "XO2", "XO3", "X2O", "X2O3", "X2O5"]
    else:
        raise ValueError("Unknown setname")

    progress_bar = tqdm.tqdm(wf_nodes)
    for node in progress_bar:
        measure_precision_out = node.outputs.measure.precision
        element = node.extras['element']

        for conf in conf_list:
            # Extract volumes and energies for each configuration
            try:
                volume_energy_dict = measure_precision_out[conf].eos.output_volume_energy.get_dict()
            except:
                print(f"WARNING! Missing output for {element} {conf}.")
                # If the configuration was not computed
                continue
            volumes = volume_energy_dict["volumes"]
            energies = volume_energy_dict["energies"]

            num_atoms = volume_energy_dict["num_of_atoms"]

            # sort
            energies = [e for _, e in sorted(zip(volumes, energies))]
            volumes = sorted(volumes)
            # List as I need to JSON-serialize it
            eos_data = (np.array([volumes, energies]).T).tolist()

            # I need to pass a numpy array
            try:
                min_volume, E0, bulk_modulus_internal, bulk_deriv, residuals = BM(np.array(eos_data))
                bulk_modulus_GPa = bulk_modulus_internal * echarge * 1.0e21
                #1 eV/Angstrom3 = 160.21766208 GPa
                bulk_modulus_ev_ang3 = bulk_modulus_GPa / 160.21766208
                BM_fit_data = {
                    'min_volume': min_volume,
                    'E0': E0,
                    'bulk_modulus_ev_ang3': bulk_modulus_ev_ang3,
                    'bulk_deriv': bulk_deriv,
                    'residuals': residuals[0]
                }
                if residuals[0] > 1.e-3:
                    warning_lines.append(f"WARNING! High fit residuals: {residuals[0]} for {element} {conf}.")
            except ValueError:
                # If we cannot find a minimum
                # Note that BM_fit_data was already set to None at the top
                warning_lines.append(f"WARNING! Unable to fit for {element} {conf}.")


            if args.setname == "unaries":
                configuration = f"X/{conf}"
            else:
                configuration = conf
            
            all_eos_data[f'{element}-{configuration}'] = eos_data
            num_atoms_in_sim_cell[f'{element}-{configuration}'] = num_atoms
            all_BM_fit_data[f'{element}-{configuration}'] = BM_fit_data
            

    data = {
        'script_version': __version__,
        # Mapping from strings like "He-X2O" to a dictionary with the UUIDs of the structure and the EOS workflow
        'uuid_mapping': uuid_mapping,
        # A list of dictionaries with information on the workchains that did not finish with a 0 exit code
        'failed_wfs': failed_wfs,
        # A dictionary that indicate for which elements and configurations there are missing outputs,
        # (only for the workchains that still had enough volumes to be considered for a fit)
        'missing_outputs': all_missing_outputs,
        # A list of dictionaries that indicate which elements and configurations have been computed completely
        # off-centre (meaning that the minimum of all computed energies is on either of the two edges, i.e. for
        # the smallest or largest volume)
        'completely_off': completely_off,
        # Dictionary with the EOS data (volumes and energies datapoints). The keys are the same as the `uuid_mapping`.
        # Values can be None.
        'eos_data': all_eos_data,
        'stress_data': all_stress_data,
        # Birch-Murnaghan fit data. See above for the keys. Can be None.
        'BM_fit_data': all_BM_fit_data,
        'num_atoms_in_sim_cell': num_atoms_in_sim_cell
    }

    # Output results to file
    os.makedirs("outputs", exist_ok=True)
    fname = f"outputs/results-{args.setname}-verification-PBE-v1-quantum_espresso-{args.output}.json"
    with open(fname, 'w') as fhandle:
        json.dump(data, fhandle, indent=2, sort_keys=True)
        
    print(f"Output written to {fname}")
    

if __name__ == "__main__":
    main()