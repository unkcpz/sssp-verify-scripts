
import os

ELEMENTS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 
            'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 
            'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 
            'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 
            'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ce', 
            'La', 'Nd', 'Pr', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 
            'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 
            'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf']

def main():
    lib_path = ".."
    dry_run = False
    FILENAMES = [i for i in os.listdir(lib_path) if 'upf' in i]
    for element in ELEMENTS:
        for fn in FILENAMES:
            if f'_{element}_' in fn:
                new_name = f"{element}.nc.oncvpsp4.spms.v1.upf"
                new_path = os.path.join(lib_path, new_name)
                if dry_run:
                    print(f'{fn} -> {new_path}')
                else:
                    os.rename(os.path.join(lib_path, fn), new_path)
        

if __name__ == "__main__":
    main()