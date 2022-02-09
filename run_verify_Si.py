#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Si'

PSEUDOS_DICT = {
        'Si_ONCV_PBE-1.2.upf': {
            'dual': 4,
            'label': 'si/sg15/z=4/nc/v1.2'
        },
        # 'Si_ONCV_PBE-1.1.upf': {
        #     'dual': 4,
        #     'label': 'si/sg15/z=4/nc/v1.1'
        # },
        # 'Si_ONCV_PBE-1.0.upf': {
        #     'dual': 4,
        #     'label': 'si/sg15/z=4/nc/v1.0'
        # },
        # 'si_pbe_v1.uspp.F.UPF': {
        #     'dual': 8,
        #     'label': 'si/gbrv/z=4/us/v1'
        # },
        # 'Si.pbe-n-kjpaw_psl.0.1.UPF': {
        #     'dual': 8,
        #     'label': 'si/psl/z=4/paw/v0.1'
        # },
        # 'Si.pbe-n-kjpaw_psl.1.0.0.UPF': {
        #     'dual': 8,
        #     'label': 'si/psl(low)/z=4/paw/v1.0.0'
        # },
        # 'Si.pbe-nl-rrkjus_psl.1.0.0.UPF': {
        #     'dual': 8,
        #     'label': 'si/psl(low)/z=4/us/v1.0.0'
        # },
        # 'Si.pbe-n-rrkjus_psl.0.1.UPF': {
        #     'dual': 8,
        #     'label': 'si/psl/z=4/us/v0.1'
        # },
    }

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
