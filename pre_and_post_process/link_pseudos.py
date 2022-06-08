"""create links to pseudopotentials to folder catogriesed by elements"""

import os
import sys

def main():
    
    element = sys.argv[1]
    func = sys.argv[2]  # functional
    
    func_path = os.path.join(f'../libraries-{func}')
    dst_path = os.path.join(f'../_sssp_{func}', element)
    
    # create element folder if not exist
    if not os.path.exists(dst_path):
        os.mkdir(dst_path)
    
    for folder in [
        "NC-DOJOv4-standard",
        "NC-DOJOv4-stringent",
        "NC-SG15-ONCVPSP4",
        "PAW-JTH1.1-standard",
        "PAW-JTH1.1-stringent",
        "PAW-PSL0.x",
        "PAW-PSL1.0.0-high",
        "PAW-PSL1.0.0-low",
        "US-GBRV-1.x",
        "US-PSL0.x",
        "US-PSL1.0.0-high",
        "US-PSL1.0.0-low",
        "UNCATOGRIZED",
        "PAW-RE-Wentzcovitch",
    ]:
        psp_folder = os.path.join(func_path, folder)
        for psp_filename in os.listdir(psp_folder):
            if psp_filename.startswith(f'{element}.'):
                src = os.path.join('..', '..', psp_folder, psp_filename)
                dst = os.path.join(dst_path, psp_filename)
                print(f'Linking: {src} -> {dst}')
                
                try:
                    os.symlink(src, dst)
                except FileExistsError:
                    os.remove(dst)
                    os.symlink(src, dst)
                    continue
                
if __name__ == '__main__':
    main()