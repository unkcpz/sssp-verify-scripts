# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "PAW-PSL-0.x"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-PSL0.x'


if mode == "pseudos_generate":
    pd_elements_lst = [
        'H', 'He', 
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se','Br', 'Kr',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd','Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
        'Cs', 'Ba',      'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt','Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra',
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb','Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf','Es', 'Fm', 'Md', 'No', 'Lr',
    ]

    def get_pseudo_list(path: Path):
        pseudo_list = []
        for element in pd_elements_lst:
            for pseudo_path in path.iterdir():
                if pseudo_path.is_file() and pseudo_path.name.split('.')[0] == element:
                    pseudo_list.append(pseudo_path.name)

        return pseudo_list

    pseudos = get_pseudo_list(Path(lib_path))

    print(len(pseudos))
    pprint.pprint(pseudos)

else:
    pseudos = [
 'H.paw.z_1.ld1.psl.v0.1.upf',
 'Li.paw.z_3.ld1.psl.v0.2.1.upf',
 'B.paw.z_3.ld1.psl.v0.1.upf',
 'C.paw.z_4.ld1.psl.v0.1.upf',
 'N.paw.z_5.ld1.psl.v0.1.upf',
 'O.paw.z_6.ld1.psl.v0.1.upf',
 'F.paw.z_7.ld1.psl.v0.1.upf',
 'Na.paw.z_9.ld1.psl.v0.2.upf',
 'Al.paw.z_3.ld1.psl.v0.1.upf',
 'Si.paw.z_4.ld1.psl.v0.1.upf',
 'P.paw.z_5.ld1.psl.v0.1.upf',
 'S.paw.z_6.ld1.psl.v0.1.upf',
 'Cl.paw.z_7.ld1.psl.v0.3.0.upf',
 'Ar.paw.z_8.ld1.psl.v0.3.0.upf',
 'Sc.paw.z_11.ld1.psl.v0.2.3.upf',
 'Fe.paw.z_16.ld1.psl.v0.2.1.upf',
 'Ni.paw.z_10.ld1.psl.v0.1.upf',
 'Cu.paw.z_11.ld1.psl.v0.2.upf',
 'Ga.paw.z_13.ld1.psl.v0.2.upf',
 'Ge.paw.z_14.ld1.psl.v0.3.1.upf',
 'As.paw.z_5.ld1.psl.v0.2.upf',
 'Se.paw.z_6.ld1.psl.v0.2.upf',
 'Br.paw.z_7.ld1.psl.v0.2.upf',
 'Sr.paw.z_10.ld1.psl.v0.2.3.upf',
 'Zr.paw.z_12.ld1.psl.v0.2.3.upf',
 'Nb.paw.z_13.ld1.psl.v0.3.0.upf',
 'Mo.paw.z_14.ld1.psl.v0.3.0.upf',
 'Tc.paw.z_15.ld1.psl.v0.3.0.upf',
 'Rh.paw.z_17.ld1.psl.v0.3.0.upf',
 'Pd.paw.z_10.ld1.psl.v0.3.0.upf',
 'Ag.paw.z_11.ld1.psl.v0.1.upf',
 'Cd.paw.z_12.ld1.psl.v0.3.1.upf',
 'In.paw.z_13.ld1.psl.v0.2.2.upf',
 'Ta.paw.z_13.ld1.psl.v0.2.upf',
 'Ir.paw.z_9.ld1.psl.v0.2.3.upf',
 'Pt.paw.z_10.ld1.psl.v0.1.upf',
 'Au.paw.z_11.ld1.psl.v0.3.0.upf',
 'Hg.paw.z_12.ld1.psl.v0.2.2.upf',
 'Tl.paw.z_13.ld1.psl.v0.2.3.upf',
 'Pb.paw.z_14.ld1.psl.v0.2.2.upf',
 'Bi.paw.z_15.ld1.psl.v0.2.2.upf',
]

    computer = 'eiger-mc-mr33-mem'
    #computer = 'daint-mc-mrcloud-mem'
    mpiprocs = 128 # 128 for eiger
    npool = 8 # 8 for eiger
    base_path = lib_path

    for pseudo in pseudos:
        pseudo_path = os.path.join(base_path, pseudo)
        command = f"aiida-sssp-workflow launch --property convergence --pw-code pw-7.0@{computer} --ph-code ph-7.0@{computer} --protocol acwf --cutoff-control standard --criteria {criteria} --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment}  -- {pseudo_path}"
        os.system(command)
        #print(command)
        print(f"Launched {pseudo}")
