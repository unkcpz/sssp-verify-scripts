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
    "Au.pbe-n-rrkjus_psl.1.0.0.UPF": {
        "dual": 8,                  
        "label": "au/psl/z=11/us/v1.0.0"     
    },                                
    "Au.pbe-n-kjpaw_psl.1.0.0.UPF": {
        "dual": 8,                  
        "label": "au/psl/z=11/paw/v1.0.0"     
    },                                  
    "Au_ONCV_PBE-1.0.upf": {       
        "dual": 4,                 
        "label": "au/sg15/z=19/nc/v1.0"
    },                      
    "Au.dojo-sr-04-std.upf": {      
        "dual": 4,                  
        "label": "au/dojo/z=19/nc/v04"   
    },                                  
    # "Au.pbe-dn-kjpaw_psl.0.3.0.UPF": {
    #     "dual": 8,                  
    #     "label": "au/psl/z=11/paw/v0.3.0"   
    # },                                   
    # "Au.pbe-spn-rrkjus_psl.1.0.0.UPF": {
    #     "dual": 8,                 
    #     "label": "au/psl(low)/z=19/us/v1.0.0"
    # },   
    #     "Au_ONCV_PBE-1.2.upf": {
    #     "dual": 4,
    #     "label": "au/sg15/z=19/nc/v1.2"
    # },
    # "Au.pbe-spfn-kjpaw_psl.1.0.0.UPF": {
    #     "dual": 8,
    #     "label": "au/psl(high)/z=33/paw/v1.0.0"
    # },
    # "Au.pbe-spfn-rrkjus_psl.1.0.0.UPF": {
    #     "dual": 8,
    #     "label": "au/psl(high)/z=33/us/v1.0.0"
    # },
    # "au_pbe_v1.uspp.F.UPF": {
    #     "dual": 8,
    #     "label": "au/gbrv/z=11/us/v1"
    # },
    # "Au.pbe-spn-kjpaw_psl.1.0.0.UPF": {
    #     "dual": 8,
    #     "label": "au/psl(low)/z=19/paw/v1.0.0"
    # },
    # "Au.pbe-dn-rrkjus_psl.0.3.0.UPF": {
    #     "dual": 8,
    #     "label": "au/psl/z=11/us/v0.3.0"
    # }
}                                                                                                            


if __name__ == '__main__':
    from aiida.orm import load_code

    pw_code = load_code('pw-6.8@eiger-hq')
    ph_code = load_code('ph-6.8@eiger-hq')

    verify_pseudos_in_folder(SSSP_DIR, ELEMENT, PSEUDOS_DICT, pw_code, ph_code, test_mode=False)
