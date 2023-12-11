# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "US-PSL-0.x"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-PSL0.x'


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
 'H.us.z_1.ld1.psl.v0.1.upf',
 'Li.us.z_3.ld1.psl.v0.2.1.upf',
 'B.us.z_3.ld1.psl.v0.1.upf',
 'C.us.z_4.ld1.psl.v0.1.upf',
 'N.us.z_5.ld1.psl.v0.1.upf',
 'O.us.z_6.ld1.psl.v0.1.upf',
 'F.us.z_7.ld1.psl.v0.1.upf',
 'Na.us.z_9.ld1.psl.v0.2.upf',
 'Al.us.z_3.ld1.psl.v0.1.upf',
 'Si.us.z_4.ld1.psl.v0.1.upf',
 'P.us.z_5.ld1.psl.v0.1.upf',
 'S.us.z_6.ld1.psl.v0.1.upf',
 'Cl.us.z_7.ld1.psl.v0.3.0.upf',
 'Ar.us.z_8.ld1.psl.v0.3.0.upf',
 'Sc.us.z_11.ld1.psl.v0.2.3.upf',
 'Fe.us.z_16.ld1.psl.v0.2.1.upf',
 'Ni.us.z_10.ld1.psl.v0.1.upf',
 'Cu.us.z_11.ld1.psl.v0.2.upf',
 'Ga.us.z_13.ld1.psl.v0.2.upf',
 'Ge.us.z_14.ld1.psl.v0.3.1.upf',
 'As.us.z_5.ld1.psl.v0.2.upf',
 'Se.us.z_6.ld1.psl.v0.2.upf',
 'Br.us.z_7.ld1.psl.v0.2.upf',
 'Sr.us.z_10.ld1.psl.v0.2.3.upf',
 'Zr.us.z_12.ld1.psl.v0.2.3.upf',
 'Nb.us.z_13.ld1.psl.v0.3.0.upf',
 'Mo.us.z_14.ld1.psl.v0.3.0.upf',
 'Tc.us.z_15.ld1.psl.v0.3.0.upf',
 'Rh.us.z_17.ld1.psl.v0.3.0.upf',
 'Pd.us.z_10.ld1.psl.v0.3.0.upf',
 'Ag.us.z_11.ld1.psl.v0.1.upf',
 'Cd.us.z_12.ld1.psl.v0.3.1.upf',
 'In.us.z_13.ld1.psl.v0.2.2.upf',
 'Ta.us.z_13.ld1.psl.v0.2.upf',
 'Ir.us.z_9.ld1.psl.v0.2.3.upf',
 'Pt.us.z_10.ld1.psl.v0.1.upf',
 'Au.us.z_11.ld1.psl.v0.3.0.upf',
 'Hg.us.z_12.ld1.psl.v0.2.2.upf',
 'Tl.us.z_13.ld1.psl.v0.2.3.upf',
 'Pb.us.z_14.ld1.psl.v0.2.2.upf',
 'Bi.us.z_15.ld1.psl.v0.2.2.upf',
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
