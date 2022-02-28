#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Ba'

PSEUDOS_DICT = {
    "Ba-sp.upf": {
        "label": "ba/dojo(optimize)-f/z=10/nc/vx"
    },
    "Ba.pbe-spn-kjpaw_psl.1.0.0.UPF": {
        "label": "ba/psl/z=10/paw/v1.0.0"
    }
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
