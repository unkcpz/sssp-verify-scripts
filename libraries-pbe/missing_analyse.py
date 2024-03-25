import sys
import os
import hashlib

# Complete list of elements up to element 118
all_elements = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", 
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", 
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", 
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", 
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", 
    "Pa", "U", "Np", "Pu", "Am", "Cm", 
]

def statistic_lib(path):
    # read the arg from the command line, the folder path list the files in the folder
    files = os.listdir(path)
    
    have_elements = []
    # for each file in the folder
    for file in files:
        # if extension is not .upf
        if not file.endswith(".upf"):
            continue

        have_elements.append(file.split(".")[0])
        
    missing_elements = list(set(all_elements) - set(have_elements))

    # find the duplicates in have_elements
    duplicates = set([x for x in have_elements if have_elements.count(x) > 1])

    if duplicates:
        print("There are duplicates in the folder")
        print(duplicates)

    print(f"In total there are {len(have_elements)} elements in the folder")

    print("The missing elements are:")
    print(', '.join(sorted(missing_elements)))

    print()
    print(f"The included elements are:")
    print(', '.join(sorted(have_elements)))

    return have_elements

def compare_two_libs(path1, path2, recollect_folder=None):
    path1_have_elements = statistic_lib(path1)
    path2_have_elements = statistic_lib(path2)

    # if path2 have elements that path1 does not have, print them out
    if set(path2_have_elements) - set(path1_have_elements):
        print("The missing elements in the first library are:")
        print(list(set(path2_have_elements) - set(path1_have_elements)))
    
    pp_new = []
    for file1 in os.listdir(path1):
        for file2 in os.listdir(path2):
            ele1 = file1.split(".")[0]
            ele2 = file2.split(".")[0]
            if ele1 != ele2:
                continue
            
            # check the md5sum of the file if they are the same element
            # Get the md5sum of the file using python API
            md5sum1 = hashlib.md5(open(os.path.join(path1, file1), 'rb').read()).hexdigest()
            md5sum2 = hashlib.md5(open(os.path.join(path2, file2), 'rb').read()).hexdigest()            
            if md5sum1 != md5sum2:
                print(f"{file2} is new")
                pp_new.append(file2)

    print(f"There are {len(pp_new)} new elements in the second library")

    if recollect_folder:
        for file in pp_new:
            os.system(f"cp {os.path.join(path2, file)} {os.path.join(recollect_folder, file)}")
            print(f"create {file} in {recollect_folder}")

def run():
    if len(sys.argv) == 2:
        statistic_lib(sys.argv[1])
    elif len(sys.argv) > 2:
        try:
            recollect_folder = sys.argv[3]
        except:
            recollect_folder = None
        else:
            if not os.path.exists(recollect_folder):
                os.makedirs(recollect_folder)
            
        compare_two_libs(sys.argv[1], sys.argv[2], recollect_folder)

if __name__ == "__main__":
    run()
