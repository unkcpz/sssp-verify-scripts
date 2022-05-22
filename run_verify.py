#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os
import click

from aiida import orm

from aiida_sssp_workflow.workflows.verifications import DEFAULT_PROPERTIES_LIST, DEFAULT_CONVERGENCE_PROPERTIES_LIST

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

if __name__ == '__main__':
    from aiida.orm import load_code
    
    psp_path = sys.argv[1]
    
    clean_workdir_level = orm.Int(1)
        
    if mode == 'test':
        pw_code = load_code('pw-6.7@localhost')
        ph_code = load_code('ph-6.7@localhost')
        protocol = orm.Str('test')
        cutoff_control = orm.Str('test')
        criteria = orm.Str('efficiency')
        option = orm.Dict(
            dict={
                "resources": {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": 1,
                },
                "max_wallclock_seconds": 1800,
                "withmpi": True,
            }
        )
        parollization = orm.Dict(dict={})
        properties_list = DEFAULT_PROPERTIES_LIST
    
    if mode == 'precheck':
        pw_code = load_code('pw-7.0@eiger-mc-mr0')
        ph_code = load_code('ph-7.0@eiger-mc-mr0')
        protocol = orm.Str('acwf')
        cutoff_control = orm.Str('precheck')
        criteria = orm.Str('efficiency')
        option = orm.Dict(
            dict={
                "resources": {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": 128,
                },
                "max_wallclock_seconds": 1800,
                "withmpi": True,
            }
        )
        parollization = orm.Dict(dict={'npool': 16})
        properties_list = DEFAULT_CONVERGENCE_PROPERTIES_LIST
        
    if mode == 'standard':
        pw_code = load_code('pw-7.0@eiger-mc-mr0')
        ph_code = load_code('ph-7.0@eiger-mc-mr0')
        protocol = orm.Str('acwf')
        cutoff_control = orm.Str('standard')
        criteria = orm.Str('efficiency')
        option = orm.Dict(
            dict={
                "resources": {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": 128,
                },
                "max_wallclock_seconds": 1800,
                "withmpi": True,
            }
        )
        parollization = orm.Dict(dict={'npool': 16})
        properties_list = DEFAULT_PROPERTIES_LIST
        


    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code, test_mode=False)
