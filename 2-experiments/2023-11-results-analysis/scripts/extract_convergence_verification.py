"""Extract JSON data from convergence verification calculations."""
import os
import json
import tqdm
from aiida import orm
from aiida.common import exceptions

from argparse import ArgumentParser

__version__ = '0.1.0'

def main():
    parser = ArgumentParser()
    parser.add_argument("libname", help="Name of the pseudo library of convergence verification calculations.")
    parser.add_argument("criteria", help="Name of the criteria to extract, e.g. 'precision' or 'efficiency'")
    parser.add_argument("output", help="Name of the output JSON file.")
    args = parser.parse_args()
    
    WORKFLOWS_GROUP_LABEL = f'{args.libname}/convergence/{args.criteria}'

    # Get all nodes in the output group
    group_node_query = orm.QueryBuilder().append(
        orm.Group, filters={'label': WORKFLOWS_GROUP_LABEL}, tag='groups',
    ).append(orm.Node, project='*', with_group='groups')
    group_node_query.distinct()
    wf_nodes = group_node_query.all(flat=True)
    
    # Extract data from each workflow, {"<element>": {..}}
    data = {}

    progress_bar = tqdm.tqdm(wf_nodes)
    for node in progress_bar:
        filename = node.inputs.pseudo.filename
        md5 = node.inputs.pseudo.md5
        element = node.extras['element']

        try:
            convergence_out = node.outputs.convergence
            max_cutoff_wfc = 0
            max_cutoff_rho = 0
            for value in convergence_out.values():
                cutoff_wfc = value.output_parameters['wavefunction_cutoff']
                cutoff_rho = value.output_parameters['chargedensity_cutoff']

                max_cutoff_wfc = max(max_cutoff_wfc, cutoff_wfc)
                max_cutoff_rho = max(max_cutoff_rho, cutoff_rho)

            element_data = {
                "filename": filename,
                "md5": md5,
                "cutoff_wfc": max_cutoff_wfc,
                "cutoff_rho": max_cutoff_rho,
            }
        except exceptions.NotExistentAttributeError:
            element_data = {
                "filename": filename,
                "md5": md5,
            }


        data[element] = element_data

    os.makedirs("outputs", exist_ok=True)
    with open(f"outputs/{args.output}.json", 'w') as f:
        json.dump(data, f, indent=4, sort_keys=True)

    
if __name__ == "__main__":
    main()