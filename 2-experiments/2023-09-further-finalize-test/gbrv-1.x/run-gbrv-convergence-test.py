# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "GBRV-1.x"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-GBRV-1.x'


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
# 2023-12-01: I did run for Au/Li/Mg/Ca/As/Mn/Ge
    pseudos = [
 #'H.us.z_1.uspp.gbrv.v1.4.upf',
 #'Li.us.z_3.uspp.gbrv.v1.4.upf',
 #'Be.us.z_4.uspp.gbrv.v1.4.upf',
 #'B.us.z_3.uspp.gbrv.v1.4.upf',
 #'C.us.z_4.uspp.gbrv.v1.2.upf',
 #'N.us.z_5.uspp.gbrv.v1.2.upf',
 #'O.us.z_6.uspp.gbrv.v1.2.upf',
 #'F.us.z_7.uspp.gbrv.v1.4.upf',
 #'Na.us.z_9.uspp.gbrv.v1.5.upf',
 #'Mg.us.z_10.uspp.gbrv.v1.4.upf',
 #'Al.us.z_3.uspp.gbrv.v1.upf',
 #'Si.us.z_4.uspp.gbrv.v1.upf',
 #'P.us.z_5.uspp.gbrv.v1.5.upf',
 #'S.us.z_6.uspp.gbrv.v1.4.upf',
 #'Cl.us.z_7.uspp.gbrv.v1.4.upf',
 #'K.us.z_9.uspp.gbrv.v1.4.upf',
 #'Ca.us.z_10.uspp.gbrv.v1.upf',
 #'Sc.us.z_11.uspp.gbrv.v1.upf',
 #'Ti.us.z_12.uspp.gbrv.v1.4.upf',
 #'V.us.z_13.uspp.gbrv.v1.4.upf',
 #'Cr.us.z_14.uspp.gbrv.v1.5.upf',
 #'Mn.us.z_15.uspp.gbrv.v1.5.upf',
 #'Fe.us.z_16.uspp.gbrv.v1.5.upf',
 #'Co.us.z_17.uspp.gbrv.v1.2.upf',
 #'Ni.us.z_18.uspp.gbrv.v1.4.upf',
 #'Cu.us.z_19.uspp.gbrv.v1.2.upf',
 #'Zn.us.z_20.uspp.gbrv.v1.upf',
 #'Ga.us.z_19.uspp.gbrv.v1.4.upf',
 #'Ge.us.z_14.uspp.gbrv.v1.4.upf',
 #'As.us.z_5.uspp.gbrv.v1.upf',
 #'Se.us.z_6.uspp.gbrv.v1.upf',
 #'Br.us.z_7.uspp.gbrv.v1.4.upf',
 #'Rb.us.z_9.uspp.gbrv.v1.upf',
 #'Sr.us.z_10.uspp.gbrv.v1.upf',
 #'Y.us.z_11.uspp.gbrv.v1.4.upf',
 #'Zr.us.z_12.uspp.gbrv.v1.upf',
 #'Nb.us.z_13.uspp.gbrv.v1.upf',
 #'Mo.us.z_14.uspp.gbrv.v1.upf',
 #'Tc.us.z_15.uspp.gbrv.v1.upf',
 #'Ru.us.z_16.uspp.gbrv.v1.2.upf',
 #'Rh.us.z_15.uspp.gbrv.v1.4.upf',
 #'Pd.us.z_16.uspp.gbrv.v1.4.upf',
 #'Ag.us.z_19.uspp.gbrv.v1.4.upf',
 #'Cd.us.z_12.uspp.gbrv.v1.upf',
 #'In.us.z_13.uspp.gbrv.v1.4.upf',
 #'Sn.us.z_14.uspp.gbrv.v1.4.upf',
 #'Sb.us.z_15.uspp.gbrv.v1.4.upf',
 #'Te.us.z_6.uspp.gbrv.v1.upf',
 #'I.us.z_7.uspp.gbrv.v1.upf',
 #'Cs.us.z_9.uspp.gbrv.v1.upf',
 #'Ba.us.z_10.uspp.gbrv.v1.upf',
 #'Hf.us.z_12.uspp.gbrv.plus4_v1.upf',
 #'Hf.us.z_12.uspp.gbrv.v1.upf',
 #'Ta.us.z_13.uspp.gbrv.v1.upf',
 #'W.us.z_14.uspp.gbrv.v1.2.upf',
 #'Re.us.z_15.uspp.gbrv.v1.2.upf',
 #'Os.us.z_16.uspp.gbrv.v1.2.upf',
 #'Ir.us.z_15.uspp.gbrv.v1.2.upf',
 #'Pt.us.z_16.uspp.gbrv.v1.4.upf',
 #'Au.us.z_11.uspp.gbrv.v1.upf',
 #'Hg.us.z_12.uspp.gbrv.v1.upf',
 #'Tl.us.z_13.uspp.gbrv.v1.2.upf',
 #'Pb.us.z_14.uspp.gbrv.v1.upf',
 #'Bi.us.z_15.uspp.gbrv.v1.upf',
 #'La.us.z_11.uspp.gbrv.v1.upf',
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
