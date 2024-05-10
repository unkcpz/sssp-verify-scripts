"""Regenerate atompaw pseudos from old UPF file
"""
import re
import sys
import subprocess
import os
import shutil
from pathlib import Path

def read_upf_content(filename):
    # Read the entire content of the UPF file
    with open(filename, 'r') as file:
        file_content = file.read()

    # Use regular expression to find the content between <PP_INFO> and </PP_INFO>
    pattern = re.compile(r'<PP_INFO>(.*?)</PP_INFO>', re.DOTALL)
    match = pattern.search(file_content)

    if match:
        # Extract the content between <PP_INFO> and </PP_INFO>, without the tags
        content = match.group(1).strip()
        return content
    else:
        return "Content not found"

def clean_inp_content(content):
    # Remove first line and whitespace of each line
    _lines = content.strip().split('\n')[1:]
    lines = []
    for line in _lines:
        # This is for wentzcovitch pseudos, in PP_INFO there are empty line
        # as comments.
        if line.strip() == '':
            break
        lines.append(line)

    cleaned_lines = [line.strip() for line in lines]
    cleaned_content = '\n'.join(cleaned_lines)

    return cleaned_content

def append_extra_content(content, extra_content):
    return content + '\n' + extra_content.strip()

def regenerate(upf_file: Path, target_folder: Path, clean=False):
    filename = upf_file.name

    if filename in [str(f) for f in target_folder.glob("*.[uU][pP][fF]")]:
        print(f"File {filename} already exists in {target_folder}")
        return
    
    content = read_upf_content(upf_file)
    content = clean_inp_content(content)

    extra_content = """
PWSCFOUT
upfdx 0.005 upfxmin -9.0 upfzmesh 1.0
0   ! END
"""
    content = append_extra_content(content, extra_content)

    # write content to file
    element = filename.split('.')[0]
    inp_filename = f"{element}.inp"

    # 2nd argument is the folder path to write the file
    inputs_folder = target_folder / "inputs"
    Path.mkdir(inputs_folder, parents=True, exist_ok=True)
    with open(inputs_folder / f"{inp_filename}", 'w') as file:
        file.write(content)

    USER = os.environ['USER']
    run_folder = target_folder / "run"
    Path.mkdir(Path(run_folder), parents=True, exist_ok=True)
    
    # copy the input file to the run folder
    os.system(f"cp {inputs_folder}/{inp_filename} {run_folder}")

    # get absolute path of run folder
    run_folder = run_folder.resolve()
    print(f"Run generation in: {str(run_folder)}")

    # run command to generate the pseudopotential
    docker_command = f"docker run -i -v {str(run_folder)}:/workdir:rw -w /workdir -u $(id -u {USER}):$(id -g {USER}) ghcr.io/pspgen/atompaw:v2024.1102 sh -c \"atompaw < {inp_filename}\""    
    print(f"Running the command: {docker_command}")

    # after running the command, copy upf file to folder and clean the run folder
    process = subprocess.Popen(docker_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # wait for the process to finish
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        print(f"Error running the command: {docker_command}")
        print(f"Error message: {stderr.decode('utf-8')}")
        sys.exit(1)
    else:
        # Find the generated file <element>.*.UPF and move it to the target folder
        # The first match
        src_path =  next(Path(run_folder).glob(f"{element}.*.UPF"))
        dst_path = target_folder / f"{filename}"

        shutil.move(src_path, dst_path)
        print(f"generate successfully")
    if clean:
        os.system(f"rm -rf {str(run_folder)}")

def main():
    if len(sys.argv) > 3 and sys.argv[3] == "clean":
        clean = True
    else:
        clean = False

    regenerate(upf_file=Path(sys.argv[1]), target_folder=Path(sys.argv[2]), clean=clean)

if __name__ == "__main__":
    main()