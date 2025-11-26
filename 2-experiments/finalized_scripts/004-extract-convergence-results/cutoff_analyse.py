#!/bin/env python

from aiida import orm
import sys
import aiida
import h5py
import pandas as pd
from tabulate import tabulate
from aiida_sssp_workflow.utils.protocol import get_protocol

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def compute_recommended_cutoffs(xs, ys, criteria) -> int:
    for x, y in zip(reversed(xs), reversed(ys)):
        cutoff = x
        if y > criteria:
            break

    return cutoff

def get_criteria(protocol, property):
    criteria = get_protocol(category='criteria', name=protocol)

    return criteria[property]['bounds'][1]


def run():
    with h5py.File('pp_verify_convergence.h5', 'r') as f:
        lib_name = 'us-gbrv-v1.x-upf2'
        protocol = 'efficiency'
        group_name = f'validate/upf/candidate/{lib_name}'

        upf_group = orm.load_group(group_name)

        properties = ['bands', 'eos', 'phonon_frequencies', 'pressure', 'cohesive_energy']
        df = pd.DataFrame(columns=properties)

        for node in upf_group.nodes:
            filename = node.filename
            print(f"dealing with {filename}")
            convergence_h5_group = f[filename]

            cutoffs = []
            for property in properties:
                print(f"on property {property}")
                criteria = get_criteria(protocol, property)

                try:
                    xs = convergence_h5_group[f'convergence_{property}']['x'][()]
                    ys = convergence_h5_group[f'convergence_{property}']['y'][()]
                except KeyError:
                    # XXX: log me to a warning file
                    print(f"not able to get {property} of {filename}")
                    cutoff = None
                else:
                    cutoff = compute_recommended_cutoffs(xs, ys, criteria)
                
                cutoffs.append(cutoff)

            df.loc[filename] = cutoffs

        table = tabulate(df, headers='keys', tablefmt='pretty')
        # XXX: table to convergence-efficiency.out and convergence-precision.out
        # the max cutoff and compare. 
        print(table)

if __name__ == '__main__':
    aiida.load_profile()

    run()
