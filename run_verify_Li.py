#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Li'

PSEUDOS_DICT = {
    "li_pbe_v1.4.uspp.F.UPF": {
        "dual": 8,
        "label": "li/gbrv/z=3/us/v1.4"
    },
    "Li.pbe-sl-kjpaw_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "li/psl(low)/z=3/paw/v1.0.0"
    },
    "Li_ONCV_PBE-1.0.upf": {
        "dual": 4,
        "label": "li/sg15/z=3/nc/v1.0"
    },
    "Li.pbe-s-kjpaw_psl.0.2.1.UPF": {
        "dual": 8,
        "label": "li/psl/z=3/paw/v0.2.1"
    },
    "Li_ONCV_PBE-1.2.upf": {
        "dual": 4,
        "label": "li/sg15/z=3/nc/v1.2"
    },
    "Li.dojo-sr-04-std.upf": {
        "dual": 4,
        "label": "li/dojo/z=3/nc/v04"
    },
    "Li.pbe-s-kjpaw_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "li/psl/z=3/paw/v1.0.0"
    },
    "Li.pbe-s-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "li/psl/z=3/us/v1.0.0"
    },
    "Li.pbe-s-rrkjus_psl.0.2.1.UPF": {
        "dual": 8,
        "label": "li/psl/z=3/us/v0.2.1"
    },
    "Li.pbe-sl-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "li/psl(low)/z=3/us/v1.0.0"
    }
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
