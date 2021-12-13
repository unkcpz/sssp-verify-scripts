#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Cs'

PSEUDOS_DICT = {
    'cs_pbe_v1.uspp.F.UPF': {
        'dual': 8,
        'label': 'cs/gbrv/z=9/us/v1'
    },
    'Cs.dojo-sr-04-std.upf': {
        'dual': 4,
        'label': 'cs/dojo(std)/z=9/nc/v04'
    },
    'Cs_ONCV_PBE-1.0.upf': {
        'dual': 4,
        'label': 'cs/sg15/z=9/nc/v1.0'
    },
    'Cs_ONCV_PBE-1.2.upf': {
        'dual': 4,
        'label': 'cs/sg15/z=9/nc/v1.2'
    },
    'Cs.pbe-spn-kjpaw_psl.1.0.0.UPF': {
        'dual': 8,
        'label': 'cs/psl/z=9/paw/v1.0.0'
    },
    'Cs.pbe-spn-kjpaw_psl.web.1.0.0.UPF': {
        'dual': 8,
        'label': 'cs/psl(web)/z=9/paw/v1.0.0'
    },
    'Cs.pbe-spnl-kjpaw_psl.1.0.0.UPF': {
        'dual': 8,
        'label': 'cs/psl(low)/z=9/paw/v1.0.0'
    },
    'Cs.pbe-spnl-rrkjus_psl.1.0.0.UPF': {
        'dual': 8,
        'label': 'cs/psl(low)/z=9/us/v1.0.0'
    },
    'Cs.pbe-spn-rrkjus_psl.1.0.0.UPF': {
        'dual': 8,
        'label': 'cs/psl/z=9/us/v1.0.0'
    },
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code, test_mode=False)
