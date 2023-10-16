"""create links to pseudopotentials to folder catogriesed by elements"""

from fnmatch import fnmatch
import os
import re
import sys
import filecmp
from turtle import st

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
        "NC-SPMS",
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
        "PAW-ACT-UNIMARBURG",
    ]:
        psp_folder = os.path.join(func_path, folder)
        for psp_filename in os.listdir(psp_folder):
            if psp_filename.startswith(f'{element}.'):
                src = os.path.join('..', '..', psp_folder, psp_filename)
                dst = os.path.join(link_dst_path, psp_filename)
                
                if folder == 'NC-DOJOv4-stringent':
                    str_src = os.path.join('.', psp_folder, psp_filename)
                    try:
                        fn = [f for f in os.listdir(os.path.join('.', func_path, "NC-DOJOv4-standard")) 
                            if fnmatch(f, f'{element}.nc.*.dojo.*.upf')][0]
                    except:
                        print(f"{fn} not exist.")
                        continue
                        
                    std_src = os.path.join('.', func_path, "NC-DOJOv4-standard", fn)
                    
                    if filecmp.cmp(str_src, std_src):
                        print(f'DOJO std and str for {element} are the same. SKIP')
                        continue
                    
                if folder == 'PAW-JTH1.1-stringent':
                    str_src = os.path.join('.', psp_folder, psp_filename)
                    
                    try:
                        fn = [f for f in os.listdir(os.path.join('.', func_path, "PAW-JTH1.1-standard")) 
                            if fnmatch(f, f'{element}.paw.*.jth.*.upf')][0]
                    except:
                        print(f"{fn} not exist.")
                        continue
                        
                    std_src = os.path.join('.', func_path, "PAW-JTH1.1-standard", fn)
                    
                    if filecmp.cmp(str_src, std_src):
                        print(f'JTH std and str for {element} are the same. SKIP')
                        continue
                
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