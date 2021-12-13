#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Fe'

PSEUDOS_DICT = {
    "Fe_ONCV_PBE-1.2.upf": {
        "dual": 4,
        "label": "fe/sg15/z=16/nc/v1.2"
    },
    "Fe.dojo-sr-04-std.upf": {
        "dual": 4,
        "label": "fe/dojo/z=16/nc/v04"
    },
    "Fe.pbe-n-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "fe/psl(low)/z=8/us/v1.0.0"
    },
    "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF": {
        "dual": 8,
        "label": "fe/psl/z=16/paw/v0.2.1"
    },
    "fe_pbe_v1.5.uspp.F.UPF": {
        "dual": 8,
        "label": "fe/gbrv/z=16/us/v1.5"
    },
    "Fe.pbe-n-kjpaw_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "fe/psl(low)/z=8/paw/v1.0.0"
    },
    "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF": {
        "dual": 8,
        "label": "fe/psl/z=16/us/v0.2.1"
    },
    "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "fe/psl/z=16/us/v1.0.0"
    },
    "Fe_ONCV_PBE-1.0.upf": {
        "dual": 4,
        "label": "fe/sg15/z=16/nc/v1.0"
    },
    "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "fe/psl/z=16/paw/v1.0.0"
    }
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
