"""Get the pseudos_dict from folder"""
import json
import os

from pseudo_parser.upf_parser import parse

def get_pseudos_dict(element_dir):
    d = {}
    
    for filename in os.listdir(element_dir):
        abspath = os.path.join(element_dir, filename)
        with open(abspath, 'r') as h:
            content = h.read()
            upf_info = parse(content)
            
        dual = 4
        pp_family = 'unknown'
        pp_type = upf_info['pp_type']
        z = upf_info['z_valence']
        element = upf_info['element']
        if 'ONCV' in filename:
            dual = 4
            pp_family = 'sg15'
        elif 'dojo' in filename:
            dual = 4
            pp_family = 'dojo'
        elif 'psl' in filename:
            dual = 8
            pp_family = 'psl'
        elif 'uspp.F' in filename:
            dual = 8
            pp_family = 'gbrv'
        
        pp_type = pp_type.lower()
        if pp_type == 'uspp':
            pp_type = 'us'
            
        d[filename] = {
            'dual': dual,
            'label': f'{element.lower()}/{pp_family}/z={z}/{pp_type}/v',
        }
    
    return d

if __name__ == '__main__':
    import sys
    
    element = sys.argv[1]
    
    folder = f'./_sssp/{element}/'
    d = get_pseudos_dict(folder)
    print(json.dumps(d, sort_keys=False, indent=4))