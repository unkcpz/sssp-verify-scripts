# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "US-PSL-1.0.0-low"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-PSL1.0.0-low'


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
 'Li.us.z_3.ld1.psl.v1.0.0-low.upf',
 'Be.us.z_2.ld1.psl.v1.0.0-low.n.upf',
 'Be.us.z_4.ld1.psl.v1.0.0-low.sl.upf',
 'O.us.z_6.ld1.psl.v1.0.0-low.upf',
 'Na.us.z_9.ld1.psl.v1.0.0-low.upf',
 'Mg.us.z_10.ld1.psl.v1.0.0-low.upf',
 'Al.us.z_3.ld1.psl.v1.0.0-low.upf',
 'Si.us.z_4.ld1.psl.v1.0.0-low.upf',
 'P.us.z_5.ld1.psl.v1.0.0-low.upf',
 'S.us.z_6.ld1.psl.v1.0.0-low.upf',
 'Cl.us.z_7.ld1.psl.v1.0.0-low.upf',
 'Ar.us.z_8.ld1.psl.v1.0.0-low.upf',
 'V.us.z_13.ld1.psl.v1.0.0-low.upf',
 'Fe.us.z_8.ld1.psl.v1.0.0-low.upf',
 'Co.us.z_9.ld1.psl.v1.0.0-low.upf',
 'Ni.us.z_10.ld1.psl.v1.0.0-low.upf',
 'Cu.us.z_11.ld1.psl.v1.0.0-low.upf',
 'Zn.us.z_12.ld1.psl.v1.0.0-low.upf',
 'Ga.us.z_13.ld1.psl.v1.0.0-low.upf',
 'Ge.us.z_4.ld1.psl.v1.0.0-low.upf',
 'As.us.z_5.ld1.psl.v1.0.0-low.upf',
 'Se.us.z_6.ld1.psl.v1.0.0-low.upf',
 'Br.us.z_7.ld1.psl.v1.0.0-low.upf',
 'Pd.us.z_10.ld1.psl.v1.0.0-low.upf',
 'Ag.us.z_11.ld1.psl.v1.0.0-low.upf',
 'Cd.us.z_12.ld1.psl.v1.0.0-low.upf',
 'Sb.us.z_5.ld1.psl.v1.0.0-low.upf',
 'Te.us.z_6.ld1.psl.v1.0.0-low.upf',
 'I.us.z_7.ld1.psl.v1.0.0-low.upf',
 'Cs.us.z_9.ld1.psl.v1.0.0-low.upf',
 'Hf.us.z_12.ld1.psl.v1.0.0-low.upf',
 'Ta.us.z_13.ld1.psl.v1.0.0-low.upf',
 'W.us.z_14.ld1.psl.v1.0.0-low.upf',
 'Re.us.z_15.ld1.psl.v1.0.0-low.upf',
 'Os.us.z_16.ld1.psl.v1.0.0-low.upf',
 'Ir.us.z_9.ld1.psl.v1.0.0-low.n.upf',
 'Ir.us.z_17.ld1.psl.v1.0.0-low.spn.upf',
 'Pt.us.z_18.ld1.psl.v1.0.0-low.spn.upf',
 'Pt.us.z_10.ld1.psl.v1.0.0-low.n.upf',
 'Au.us.z_19.ld1.psl.v1.0.0-low.spn.upf',
 'Au.us.z_11.ld1.psl.v1.0.0-low.n.upf',
 'Hg.us.z_12.ld1.psl.v1.0.0-low.upf',
 #'Ce.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Pr.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Nd.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Pm.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Sm.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Eu.us.z_10.ld1.psl.v1.0.0-low.spn.upf',
 #'Eu.us.z_11.ld1.psl.v1.0.0-low.spdn.upf',
 #'Gd.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Tb.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Dy.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Ho.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Er.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Tm.us.z_11.ld1.psl.v1.0.0-low.upf',
 #'Yb.us.z_11.ld1.psl.v1.0.0-low.spdn.upf',
 #'Yb.us.z_10.ld1.psl.v1.0.0-low.spn.upf',
 #'Lu.us.z_11.ld1.psl.v1.0.0-low.upf',
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
