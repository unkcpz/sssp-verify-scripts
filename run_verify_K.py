#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'K'

PSEUDOS_DICT = {                             
    "K_ONCV_PBE-1.0.upf": {
        "dual": 4,                
        "label": "k/sg15/z=9/nc/v1.0"
    },                                 
    "K_ONCV_PBE-1.2.upf": {
        "dual": 4,               
        "label": "k/sg15/z=9/nc/v1.0"
    },                       
    "K.pbe-spn-kjpaw_psl.1.0.0.UPF": {
        "dual": 8,                
        "label": "k/psl/z=9/paw/v1.0.0"
    },
    "k_pbe_v1.4.uspp.F.UPF": {                                                                                         
        "dual": 8,
        "label": "k/gbrv/z=9/us/v1.4"
    },
    "K.pbe-spn-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "k/psl/z=9/us/v1.0.0"
    },
    "K.dojo-sr-04-std.upf": {
        "dual": 4,
        "label": "k/dojo/z=9/nc/v04"
    }
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
