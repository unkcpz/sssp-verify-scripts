"""Write z_valence to pseudo files of given library folder"""

import os
import sys
from pathlib import Path

from aiida.plugins import DataFactory

UpfData = DataFactory('pseudo.upf')

def main(lib_path, dry_run=True):
    print(lib_path)
    
    for filename in os.listdir(lib_path):
        psp_path = os.path.join(lib_path, filename)
        path = Path(psp_path)
        if not path.is_file():
            print(f'{path} is a folder, skip.')
            continue
        with open(psp_path, "rb") as stream:
            try:
                pseudo = UpfData(stream)
                z_valence = pseudo.z_valence
            except Exception as exc:
                print(f'{filename} got error {exc}')
        
        if '.z_' in filename:
            # raise ValueError(f'Seems already has z_valence in filename {filename}')
            print(f'Seems already has z_valence set for file {filename}')
            continue
        
        if '.nc.' in filename:
            new_name = filename.replace('.nc.', f'.nc.z_{z_valence}.')
        elif '.us.' in filename:
            new_name = filename.replace('.us.', f'.us.z_{z_valence}.')
        elif '.paw.' in filename:
            new_name = filename.replace('.paw.', f'.paw.z_{z_valence}.')
        else:
            raise ValueError(f'Not a proper pseudopotential filename {filename}')
        
        new_path = os.path.join(lib_path, new_name)
        if dry_run:
            print(f'{psp_path} -> {new_path}')
        else:
            os.rename(psp_path, new_path)

def parse(string):
    d = {'True': True, 'False': False}
    return d.get(string, string)

if __name__ == '__main__':
    import aiida
    
    aiida.load_profile()
    
    lib_path = sys.argv[1]
    
    # python psp_fn_z.py <lib> dry_run True
    dry_run = parse(sys.argv[3])
    main(lib_path, dry_run)