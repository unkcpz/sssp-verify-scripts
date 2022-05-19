# Extract input from upf file to dat file
# input are lines between `<PP_INPUTFILE>`  `</PP_INPUTFILE>`
import glob
import os

def main():
    src_dir = './input_generate'
    dst_dir = './_build'
    # clean all .dat files
    for f in glob.glob(f"{dst_dir}/*.dat"):
        os.remove(f)
        
    for f in glob.glob(f"{src_dir}/*-1.2.upf"):
        basename = os.path.basename(f)
        input_filename = basename.replace('upf', 'dat')
        with open(f'./{f}', 'r') as fh:
            lines = fh.readlines()
            for ln, line in enumerate(lines):
                if '<PP_INPUTFILE>' in line:
                    ln_start = ln
                    
                if '</PP_INPUTFILE>' in line:
                    ln_end = ln
                    break
                
            input_lines = lines[ln_start+1:ln_end]
            
            input = ''.join(input_lines)
            
        with open(os.path.join(dst_dir, input_filename), 'w') as fh:
            fh.write(input)

if __name__ == '__main__':
    main()