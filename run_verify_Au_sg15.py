#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Running verification workchain
"""
import os

from sssp_verify_scripts import verify_pseudos_in_folder

SSSP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_sssp')

ELEMENT = 'Au'

PSEUDOS_DICT = {                                                         
    "Au_ONCV_PBE-1.0.oncvpsp2.upf": {                     
        "label": "au/sg15(onvc2)/z=19/nc/v1.0"
    },       
    "Au_ONCV_PBE-1.0.oncvpsp3.upf": {                    
        "label": "au/sg15(onvc3)/z=19/nc/v1.0"
    },       
    "Au_ONCV_PBE-1.0.oncvpsp4.upf": {                    
        "label": "au/sg15(onvc4)/z=19/nc/v1.0"
    },           
    "au_pbe_v1.uspp.F.UPF": {
        "label": "au/gbrv/z=11/us/v1"
    },
    "Au.dojo-sr-04-std.upf": {      
        "dual": 4,                  
        "label": "au/dojo/z=19/nc/v04"   
    },
    "Au.pbe-spn-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,                 
        "label": "au/psl(low)/z=19/us/v1.0.0"
    },   
    "Au.pbe-spfn-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,
        "label": "au/psl(high)/z=33/us/v1.0.0"
    },
}                                                                                                            


if __name__ == '__main__':
    from aiida.orm import load_code

    # pw_code = load_code('pw-6.8@eiger-hq')
    # ph_code = load_code('ph-6.8@eiger-hq')
    pw_code = load_code('pw-6.8@eiger-mc-mr0')
    ph_code = load_code('ph-6.8@eiger-mc-mr0')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code)
