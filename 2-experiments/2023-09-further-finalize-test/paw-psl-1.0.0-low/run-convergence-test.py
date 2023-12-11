# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "PAW-PSL-1.0.0-low"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-PSL1.0.0-low'


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
 'Li.paw.z_3.ld1.psl.v1.0.0-low.upf',
 'Be.paw.z_4.ld1.psl.v1.0.0-low.sl.upf',
 'Be.paw.z_2.ld1.psl.v1.0.0-low.upf',
 'O.paw.z_6.ld1.psl.v1.0.0-low.upf',
 'Na.paw.z_9.ld1.psl.v1.0.0-low.upf',
 'Mg.paw.z_10.ld1.psl.v1.0.0-low.upf',
 'Al.paw.z_3.ld1.psl.v1.0.0-low.upf',
 'Si.paw.z_4.ld1.psl.v1.0.0-low.upf',
 'P.paw.z_5.ld1.psl.v1.0.0-low.upf',
 'S.paw.z_6.ld1.psl.v1.0.0-low.upf',
 'Cl.paw.z_7.ld1.psl.v1.0.0-low.upf',
 'Ar.paw.z_8.ld1.psl.v1.0.0-low.upf',
 'V.paw.z_13.ld1.psl.v1.0.0-low.upf',
 'Fe.paw.z_8.ld1.psl.v1.0.0-low.upf',
 'Co.paw.z_9.ld1.psl.v1.0.0-low.upf',
 'Ni.paw.z_10.ld1.psl.v1.0.0-low.upf',
 'Cu.paw.z_11.ld1.psl.v1.0.0-low.upf',
 'Zn.paw.z_12.ld1.psl.v1.0.0-low.upf',
 'Ga.paw.z_13.ld1.psl.v1.0.0-low.upf',
 'Ge.paw.z_4.ld1.psl.v1.0.0-low.upf',
 'As.paw.z_5.ld1.psl.v1.0.0-low.upf',
 'Se.paw.z_6.ld1.psl.v1.0.0-low.upf',
 'Br.paw.z_7.ld1.psl.v1.0.0-low.upf',
 'Pd.paw.z_10.ld1.psl.v1.0.0-low.upf',
 'Ag.paw.z_11.ld1.psl.v1.0.0-low.upf',
 'Cd.paw.z_12.ld1.psl.v1.0.0-low.upf',
 'Sb.paw.z_5.ld1.psl.v1.0.0-low.upf',
 'Te.paw.z_6.ld1.psl.v1.0.0-low.upf',
 'I.paw.z_7.ld1.psl.v1.0.0-low.upf',
 'Cs.paw.z_9.ld1.psl.v1.0.0-low.upf',
 'Hf.paw.z_12.ld1.psl.v1.0.0-low.upf',
 'Ta.paw.z_13.ld1.psl.v1.0.0-low.upf',
 'W.paw.z_14.ld1.psl.v1.0.0-low.upf',
 'Re.paw.z_15.ld1.psl.v1.0.0-low.upf',
 'Os.paw.z_16.ld1.psl.v1.0.0-low.upf',
 'Ir.paw.z_9.ld1.psl.v1.0.0-low.n.upf',
 'Ir.paw.z_17.ld1.psl.v1.0.0-low.spn.upf',
 'Pt.paw.z_10.ld1.psl.v1.0.0-low.n.upf',
 'Pt.paw.z_18.ld1.psl.v1.0.0-low.spn.upf',
 'Au.paw.z_11.ld1.psl.v1.0.0-low.n.upf',
 'Au.paw.z_19.ld1.psl.v1.0.0-low.spn.upf',
 'Hg.paw.z_12.ld1.psl.v1.0.0-low.upf',
 #'Ce.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Pr.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Nd.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Pm.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Sm.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Eu.paw.z_10.ld1.psl.v1.0.0-low.spn.upf',
 #'Eu.paw.z_11.ld1.psl.v1.0.0-low.spdn.upf',
 #'Gd.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Tb.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Dy.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Ho.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Er.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Tm.paw.z_11.ld1.psl.v1.0.0-low.upf',
 #'Yb.paw.z_11.ld1.psl.v1.0.0-low.spdn.upf',
 #'Yb.paw.z_10.ld1.psl.v1.0.0-low.spn.upf',
 #'Lu.paw.z_11.ld1.psl.v1.0.0-low.upf',
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
