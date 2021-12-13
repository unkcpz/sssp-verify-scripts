#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'F'

PSEUDOS_DICT = {
    "F_ONCV_PBE-1.2.upf": {
        "dual": 4,
        "label": "f/sg15/z=7/nc/v1.2"
    },
    "F.pbe-n-kjpaw_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "f/psl/z=7/paw/v1.0.0"
    },
    "F.dojo-sr-04-std.upf": {
        "dual": 4,
        "label": "f/dojo/z=7/nc/v04"
    },
    "F.pbe-n-rrkjus_psl.0.1.UPF": {
        "dual": 8,
        "label": "f/psl/z=7/us/v0.1"
    },
    "F.pbe-n-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "f/psl/z=7/us/v1.0.0"
    },
    "F_ONCV_PBE-1.0.upf": {
        "dual": 4,
        "label": "f/sg15/z=7/nc/v1.0"
    },
    "F.pbe-n-kjpaw_psl.0.1.UPF": {
        "dual": 8,
        "label": "f/psl/z=7/paw/v0.1"
    },
    "f_pbe_v1.4.uspp.F.UPF": {
        "dual": 8,
        "label": "f/gbrv/z=7/us/v1.4"
    }
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
