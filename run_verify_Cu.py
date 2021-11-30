#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Cu'

PSEUDOS_DICT = {
        'Cu_ONCV_PBE-1.2.upf': {
            'dual': 4,
            'label': 'cu/sg15/z=19/nc/v1.2'
        },
        'cu_pbe_v1.2.uspp.F.UPF': {
            'dual': 8,
            'label': 'cu/gbrv/z=19/us/v1.2'
        },
        'Cu.pbe-dn-kjpaw_psl.0.2.UPF': {
            'dual': 8,
            'label': 'cu/psl(low)/z=11/paw/v0.2'
        },
        'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF': {
            'dual': 8,
            'label': 'cu/psl(low)/z=11/paw/v1.0.0'
        },
        'Cu.pbe-dn-rrkjus_psl.0.2.UPF': {
            'dual': 8,
            'label': 'cu/psl(low)/z=11/us/v0.2'
        },
        'Cu.pbe-dn-rrkjus_psl.1.0.0.UPF': {
            'dual': 8,
            'label': 'cu/psl(low)/z=11/us/v1.0.0'
        },
        'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF': {
            'dual': 8,
            'label': 'cu/psl/z=19/pav/v1.0.0'
        }, 
        'Cu.pbe-spn-rrkjus_psl.1.0.0.UPF': {
            'dual': 8,
            'label': 'cu/psl/z=19/us/v1.0.0'
        },  
    }

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
