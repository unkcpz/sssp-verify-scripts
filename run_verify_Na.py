#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Na'

PSEUDOS_DICT = {
    "Na.pbe-spnl-rrkjus_psl.1.0.0.UPF": {
        "label": "na/psl(low)/z=9/us/v"
    },
    "Na.pbe-spn-rrkjus_psl.0.2.UPF": {
        "label": "na/psl/z=9/us/v"
    },
    "Na.dojo-sr-04-std.upf": {
        "label": "na/dojo/z=9/nc/v"
    },
    "Na.GGA-PBE-paw.UPF": {
        "label": "na/atompaw/z=9/paw/v"
    },
    "Na.pbe-spn-rrkjus_psl.1.0.0.UPF": {
        "label": "na/psl/z=9/us/v"
    },
}

if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
