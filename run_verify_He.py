#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'He'

PSEUDOS_DICT = {
    "He.dojo-sr-04-std.upf": {
        "dual": 4,
        "label": "he/dojo/z=2/nc/v04"
    },
    "He_ONCV_PBE-1.0.upf": {
        "dual": 4,
        "label": "he/sg15/z=2/nc/v1.0"
    },
    "He.pbe-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "he/psl/z=2/us/v1.0.0"
    },
    "He_ONCV_PBE-1.2.upf": {
        "dual": 4,
        "label": "he/sg15/z=2/nc/v1.2"
    },
    "He.pbe-kjpaw_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "he/psl/z=2/paw/v1.0.0"
    }
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
