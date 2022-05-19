# Extract upf content from ONCVPSP out
import enum
import glob
import os

from numpy import extract

def main():
    src_dir = './_build'
    dst_dir = '../'
    
    for f in glob.glob(f'{src_dir}/*.out'):
        basename = os.path.basename(f)
        upf_filename = basename.replace('_ONCV_PBE-1.2.dat.out', '.nc.sg15.oncvpsp4.upf')
        with open(f, 'r') as fh:
            lines = fh.readlines()
            
            for ln, line in enumerate(lines):
                if 'Begin PSP_UPF' in line:
                    ln_start = ln
                    
                if 'END_PSP' in line:
                    ln_end = ln
                    
            try:
                upf_lines = lines[ln_start+1:ln_end-1]
                upf_content = ''.join(upf_lines)
            except:
                print(f'{upf_filename} not generated.')
                continue
            
        with open(os.path.join(dst_dir, upf_filename), 'w') as fh:
            fh.write(upf_content)
 
if __name__ == '__main__':
    main()