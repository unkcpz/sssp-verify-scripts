# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "NC-LD1-v1.0.0"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-PSL1.0.0'


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
 #'H.nc.z_1.ld1.psl.v1.0.0.upf', 
 #'He.nc.z_2.ld1.psl.v1.0.0.upf',
 'Li.nc.z_1.ld1.psl.v1.0.0.upf',
 #'Be.nc.z_2.ld1.psl.v1.0.0.upf',
 #'B.nc.z_3.ld1.psl.v1.0.0.upf', 
 #'C.nc.z_4.ld1.psl.v1.0.0.upf',  
 #'N.nc.z_5.ld1.psl.v1.0.0.upf',  
 #'O.nc.z_6.ld1.psl.v1.0.0.upf',  
 #'F.nc.z_7.ld1.psl.v1.0.0.upf',  
 #'Ne.nc.z_8.ld1.psl.v1.0.0.upf', 
 'Na.nc.z_1.ld1.psl.v1.0.0.upf', 
 #'Mg.nc.z_2.ld1.psl.v1.0.0.upf', 
 #'Al.nc.z_3.ld1.psl.v1.0.0.upf', 
 #'Si.nc.z_4.ld1.psl.v1.0.0.upf', 
 #'P.nc.z_5.ld1.psl.v1.0.0.upf',  
 #'S.nc.z_6.ld1.psl.v1.0.0.upf', 
 #'Cl.nc.z_7.ld1.psl.v1.0.0.upf',
 #'Ar.nc.z_8.ld1.psl.v1.0.0.upf',
 'K.nc.z_1.ld1.psl.v1.0.0.upf', 
 'Ca.nc.z_2.ld1.psl.v1.0.0.upf',
 #'Sc.nc.z_3.ld1.psl.v1.0.0.upf',
 #'Ti.nc.z_4.ld1.psl.v1.0.0.upf', 
 #'V.nc.z_5.ld1.psl.v1.0.0.upf',  
 #'Cr.nc.z_6.ld1.psl.v1.0.0.upf',
 #'Fe.nc.z_8.ld1.psl.v1.0.0.upf',
 #'Co.nc.z_9.ld1.psl.v1.0.0.upf',
 #'Ni.nc.z_10.ld1.psl.v1.0.0.upf',
 #'Zn.nc.z_12.ld1.psl.v1.0.0.upf',
 #'Ga.nc.z_3.ld1.psl.v1.0.0.upf', 
 #'As.nc.z_5.ld1.psl.v1.0.0.upf',
 #'Se.nc.z_6.ld1.psl.v1.0.0.upf',
 #'Br.nc.z_7.ld1.psl.v1.0.0.upf',
 #'Kr.nc.z_8.ld1.psl.v1.0.0.upf',
 'Rb.nc.z_1.ld1.psl.v1.0.0.upf',
 'Sr.nc.z_2.ld1.psl.v1.0.0.upf',
 'Y.nc.z_3.ld1.psl.v1.0.0.upf', 
 #'Zr.nc.z_4.ld1.psl.v1.0.0.upf',
 #'Nb.nc.z_5.ld1.psl.v1.0.0.upf',
 #'Mo.nc.z_6.ld1.psl.v1.0.0.upf', 
 #'Tc.nc.z_7.ld1.psl.v1.0.0.upf', 
 #'Ru.nc.z_8.ld1.psl.v1.0.0.upf', 
 #'Rh.nc.z_9.ld1.psl.v1.0.0.upf', 
 #'Pd.nc.z_10.ld1.psl.v1.0.0.upf',
 #'Ag.nc.z_11.ld1.psl.v1.0.0.upf',
 #'Cd.nc.z_12.ld1.psl.v1.0.0.upf',
 #'In.nc.z_3.ld1.psl.v1.0.0.upf', 
 #'Sn.nc.z_4.ld1.psl.v1.0.0.upf', 
 #'Sb.nc.z_5.ld1.psl.v1.0.0.upf', 
 #'Te.nc.z_6.ld1.psl.v1.0.0.upf', 
 #'I.nc.z_7.ld1.psl.v1.0.0.upf',  
 'Xe.nc.z_8.ld1.psl.v1.0.0.upf', 
 'Cs.nc.z_1.ld1.psl.v1.0.0.upf',
 'Ba.nc.z_2.ld1.psl.v1.0.0.upf',
 #'Hf.nc.z_4.ld1.psl.v1.0.0.upf',
 #'Re.nc.z_7.ld1.psl.v1.0.0.upf',
 #'Os.nc.z_8.ld1.psl.v1.0.0.upf',
 #'Ir.nc.z_9.ld1.psl.v1.0.0.upf',
 #'Pt.nc.z_10.ld1.psl.v1.0.0.upf',
 #'Au.nc.z_11.ld1.psl.v1.0.0.upf',
 #'Hg.nc.z_12.ld1.psl.v1.0.0.upf',
 #'Tl.nc.z_3.ld1.psl.v1.0.0.upf',
 #'Pb.nc.z_4.ld1.psl.v1.0.0.upf',
 #'Bi.nc.z_5.ld1.psl.v1.0.0.upf',
 #'Po.nc.z_6.ld1.psl.v1.0.0.upf',
 #'At.nc.z_7.ld1.psl.v1.0.0.upf',
 #'Rn.nc.z_8.ld1.psl.v1.0.0.upf',
 'Fr.nc.z_1.ld1.psl.v1.0.0.upf',
 'Ra.nc.z_2.ld1.psl.v1.0.0.upf',
]


    #computer = 'eiger-mc-mr32-mem'
    computer = 'daint-mc-mrcloud-mem'
    mpiprocs = 36 # 128 for eiger
    npool = 4 # 8 for eiger
    base_path = lib_path

    for pseudo in pseudos:
        pseudo_path = os.path.join(base_path, pseudo)
        command = f"aiida-sssp-workflow launch --property convergence --pw-code pw-7.0@{computer} --ph-code ph-7.0@{computer} --protocol acwf --cutoff-control standard --criteria {criteria} --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment}  -- {pseudo_path}"
        os.system(command)
        #print(command)
        print(f"Launched {pseudo}")
