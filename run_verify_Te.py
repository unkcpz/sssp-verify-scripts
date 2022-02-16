#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Te'

PSEUDOS_DICT = {
    "Te.pbe-dn-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "te/psl/z=16/us/v1.0.0"
    },
    # "Te.pbe-n-rrkjus_psl.1.0.0.UPF": {
    #     "dual": 8,
    #     "label": "te/psl(low)/z=6/us/v1.0.0"
    # },
    # "Te.pbe-dn-kjpaw_psl.1.0.0.UPF": {
    #     "dual": 8,
    #     "label": "te/psl/z=16/paw/v1.0.0"
    # },
    # "Te.pbe-n-kjpaw_psl.1.0.0.UPF": {
    #     "dual": 8,
    #     "label": "te/psl(low)/z=6/paw/v1.0.0"
    # },
    # "Te_ONCV_PBE-1.0.upf": {
    #     "dual": 4,
    #     "label": "te/sg15/z=16/nc/v1.0"
    # },
    # "Te.pbe-dn-rrkjus_psl.0.3.1.UPF": {
    #     "dual": 8,
    #     "label": "te/psl/z=16/us/v0.3.1"
    # },
    "Te.dojo-sr-04-std.upf": {
        "dual": 4,
        "label": "te/dojo/z=16/nc/v04"
    },
    "Te_ONCV_PBE-1.2.upf": {
        "dual": 4,
        "label": "te/sg15/z=16/nc/v1.2"
    },
    # "Te_ONCV4_PBE-1.2.upf": {
    #     "dual": 4,
    #     "label": "te/sg15(oncv4)/z=16/nc/v1.2"
    # },
    "te_pbe_v1.uspp.F.UPF": {
        "dual": 8,
        "label": "te/gbrv/z=6/us/v1"
    },
    # "Te.pbe-dn-kjpaw_psl.0.3.1.UPF": {
    #     "dual": 8,
    #     "label": "te/psl/z=16/paw/v0.3.1"
    # }
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
