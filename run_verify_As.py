#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'As'

PSEUDOS_DICT = {
    # 'as_pbe_v1.uspp.F.UPF': {
    #     'dual': 8,
    #     'label': 'as/gbrv/z=5/us/v1'
    # },
    'As.dojo-sr-04-std.upf': {
        'dual': 4,
        'label': 'as/dojo(std)/z=15/nc/v04'
    },
    # 'As_ONCV_PBE-1.0.upf': {
    #     'dual': 4,
    #     'label': 'as/sg15/z=5/nc/v1.0'
    # },
    # 'As_ONCV_PBE-1.2.upf': {
    #     'dual': 4,
    #     'label': 'as/sg15/z=5/nc/v1.2'
    # },
    # 'As.pbe-dn-kjpaw_psl.1.0.0.UPF': {
    #     'dual': 8,
    #     'label': 'as/psl/z=15/paw/v1.0.0'
    # },
    'As.pbe-dn-rrkjus_psl.1.0.0.UPF': {
        'dual': 8,
        'label': 'as/psl/z=15/us/v1.0.0'
    },
    # 'As.pbe-n-kjpaw_psl.0.2.UPF': {
    #     'dual': 8,
    #     'label': 'as/psl(low)/z=5/paw/v0.2'
    # },
    # 'As.pbe-n-kjpaw_psl.1.0.0.UPF': {
    #     'dual': 8,
    #     'label': 'as/psl(low)/z=5/paw/v1.0.0'
    # },
    'As.pbe-n-rrkjus_psl.0.2.UPF': {
        'dual': 8,
        'label': 'as/psl(low)/z=5/us/v0.2'
    },
    # 'As.pbe-n-rrkjus_psl.1.0.0.UPF': {
    #     'dual': 8,
    #     'label': 'as/psl(low)/z=5/us/v1.0.0'
    # },
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
