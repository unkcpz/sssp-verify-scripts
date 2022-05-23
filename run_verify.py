#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os
import click

import aiida
from aiida import orm
from aiida.plugins import DataFactory

from aiida_sssp_workflow.workflows.verifications import DEFAULT_PROPERTIES_LIST, DEFAULT_CONVERGENCE_PROPERTIES_LIST

from sssp_verify_scripts import run_verification

UpfData = DataFactory('pseudo.upf')

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

def inputs_from_mode(mode):
    inputs = {}
    if mode == 'TEST':
        inputs['pw_code'] = orm.load_code('pw-6.7@localhost')
        inputs['ph_code'] = orm.load_code('ph-6.7@localhost')
        inputs['protocol'] = orm.Str('test')
        inputs['cutoff_control'] = orm.Str('test')
        inputs['criteria'] = orm.Str('efficiency')
        inputs['options'] = orm.Dict(
            dict={
                "resources": {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": 1,
                },
                "max_wallclock_seconds": 1800,
                "withmpi": False,
            }
        )
        inputs['parallization'] = orm.Dict(dict={})
        inputs['properties_list'] = DEFAULT_PROPERTIES_LIST
        
    if mode == 'PRECHECK':
        inputs['pw_code'] = orm.load_code('pw-7.0@eiger-mc-mr0')
        inputs['ph_code'] = orm.load_code('pw-7.0@eiger-mc-mr0')
        inputs['protocol'] = orm.Str('acwf')
        inputs['cutoff_control'] = orm.Str('precheck')
        inputs['criteria'] = orm.Str('efficiency')
        inputs['options'] = orm.Dict(
            dict={
                "resources": {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": 128,
                },
                "max_wallclock_seconds": 1800,
                "withmpi": True,
            }
        )
        inputs['parallization'] = orm.Dict(dict={'npool': 16})
        inputs['properties_list'] = DEFAULT_CONVERGENCE_PROPERTIES_LIST
        
    # if mode == 'standard':
    #     pw_code = load_code('pw-7.0@eiger-mc-mr0')
    #     ph_code = load_code('ph-7.0@eiger-mc-mr0')
    #     protocol = orm.Str('acwf')
    #     cutoff_control = orm.Str('standard')
    #     criteria = orm.Str('efficiency')
    #     option = orm.Dict(
    #         dict={
    #             "resources": {
    #                 "num_machines": 1,
    #                 "num_mpiprocs_per_machine": 128,
    #             },
    #             "max_wallclock_seconds": 1800,
    #             "withmpi": True,
    #         }
    #     )
    #     parollization = orm.Dict(dict={'npool': 16})
    #     properties_list = DEFAULT_PROPERTIES_LIST
        
    return inputs

@click.command()
@click.option('profile', '-p', help='profile')
@click.option('--mode', type=click.Choice(['TEST', 'PRECHECK', 'STANDARD'], case_sensitive=False), 
              help='mode of verification.')
@click.argument('filename', type=click.Path(exists=True))
def run(profile, mode, filename):
    click.echo(profile)

    aiida.load_profile(profile)
    
    inputs = inputs_from_mode(mode=mode)
    
    basename = os.path.basename(filename)
    label, _ = os.path.splitext(basename)
    label = f'({mode}) {label}'

    with open(filename, "rb") as stream:
        pseudo = UpfData(stream)
        
    node = run_verification(
        **inputs, 
        **{
            'pseudo': pseudo,
            'label': label,
        }, 
    )

    click.echo(node)
    click.echo(click.format_filename(filename))

    # verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code, test_mode=False)


if __name__ == '__main__':
    run()