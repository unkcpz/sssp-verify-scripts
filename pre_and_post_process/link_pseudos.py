"""create links to pseudopotentials to folder catogriesed by elements"""

import os
import sys

def main():
    
    element = sys.argv[1]
    func = sys.argv[2]  # functional
    
    func_path = os.path.join(f'./libraries-{func}')
    link_dst_path = os.path.join(f'./_sssp_{func}', element)
    
    # create element folder if not exist
    if not os.path.exists(os.path.join(f'./_sssp_{func}', element)):
        os.mkdir(os.path.join(f'./_sssp_{func}', element))
    
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
        "PAW-RE-Wentzcovitch/legacy",
        "PAW-RE-Wentzcovitch/neo",
    ]:
        psp_folder = os.path.join(func_path, folder)
        for psp_filename in os.listdir(psp_folder):
            if psp_filename.startswith(f'{element}.'):
                src = os.path.join('..', '..', psp_folder, psp_filename)
                dst = os.path.join(link_dst_path, psp_filename)
                
                try:
                    os.symlink(src, dst)
                except FileExistsError:
                    os.remove(dst)
                    os.symlink(src, dst)
                    print(f'Exist: Clean and Linking: {src} -> {dst}')
                    continue
                except FileNotFoundError:
                    # for example no SG15 lanthanides
                    continue
                else:
                    print(f'Linking: {src} -> {dst}')
                
if __name__ == '__main__':
    main()