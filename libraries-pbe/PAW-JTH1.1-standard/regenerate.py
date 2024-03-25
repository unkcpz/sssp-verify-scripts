import re
import sys
import subprocess
import os

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
    lines = content.strip().split('\n')[1:]
    cleaned_lines = [line.strip() for line in lines]
    cleaned_content = '\n'.join(cleaned_lines)

    return cleaned_content

def append_extra_content(content, extra_content):
    return content + '\n' + extra_content.strip()

if __name__ == "__main__":
    # CLI: get filename from command line argument
    filename = sys.argv[1]
    
    content = read_upf_content(filename)
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

    # 2nd argument is the folder to write the file
    folder = sys.argv[2]
    with open(f"{folder}/{inp_filename}", 'w') as file:
        file.write(content)

    USER = os.environ['USER']
    RUN_FOLDER = f"{folder}/run"
    if not os.path.exists(RUN_FOLDER):
        os.makedirs(RUN_FOLDER)
        # copy the input file to the run folder
    
    os.system(f"cp {folder}/{inp_filename} {RUN_FOLDER}")

    # get absolute path of run folder
    RUN_FOLDER = os.path.abspath(RUN_FOLDER)
    print(f"Run generation in: {RUN_FOLDER}")

    # run command to generate the pseudopotential
    docker_command = f"docker run -i -v {RUN_FOLDER}:/workdir:rw -w /workdir -u $(id -u {USER}):$(id -g {USER}) pspgen/atompaw:main sh -c \"atompaw < {inp_filename}\""    
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
        print(f"generate successfully, copy {element}.upf to {folder}")
        os.system(f"cp {RUN_FOLDER}/*.UPF {folder}")

        # Rename the UPF file to the original name
        os.system(f"mv {folder}/{element}.*.UPF {folder}/{filename}")
        # os.system(f"rm -rf {RUN_FOLDER}")

