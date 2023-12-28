# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "PAW-JTH1.1-standard"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-JTH1.1-standard'


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
 'H.paw.z_1.atompaw.jth.v1.1-std.upf',
 'He.paw.z_2.atompaw.jth.v1.1-std.upf',
 'Li.paw.z_3.atompaw.jth.v1.1-std.upf',
 'Be.paw.z_4.atompaw.jth.v1.1-std.upf',
 'B.paw.z_3.atompaw.jth.v1.1-std.upf', 
 'C.paw.z_4.atompaw.jth.v1.1-std.upf', 
 'N.paw.z_5.atompaw.jth.v1.1-std.upf', 
 'O.paw.z_6.atompaw.jth.v1.1-std.upf', 
 'F.paw.z_7.atompaw.jth.v1.1-std.upf', 
 'Ne.paw.z_8.atompaw.jth.v1.1-std.upf',
 'Na.paw.z_9.atompaw.jth.v1.1-std.upf', 
 'Mg.paw.z_10.atompaw.jth.v1.1-std.upf',
 'Al.paw.z_3.atompaw.jth.v1.1-std.upf', 
 'Si.paw.z_4.atompaw.jth.v1.1-std.upf', 
 'P.paw.z_5.atompaw.jth.v1.1-std.upf',  
 'S.paw.z_6.atompaw.jth.v1.1-std.upf',  
 'Cl.paw.z_7.atompaw.jth.v1.1-std.upf', 
 'Ar.paw.z_8.atompaw.jth.v1.1-std.upf', 
 'K.paw.z_9.atompaw.jth.v1.1-std.upf',  
 'Ca.paw.z_10.atompaw.jth.v1.1-std.upf',
 'Sc.paw.z_11.atompaw.jth.v1.1-std.upf',
 'Ti.paw.z_12.atompaw.jth.v1.1-std.upf',
 'V.paw.z_13.atompaw.jth.v1.1-std.upf', 
 'Cr.paw.z_14.atompaw.jth.v1.1-std.upf',
 'Mn.paw.z_15.atompaw.jth.v1.1-std.upf',
 'Fe.paw.z_16.atompaw.jth.v1.1-std.upf',
 'Co.paw.z_17.atompaw.jth.v1.1-std.upf',
 'Ni.paw.z_18.atompaw.jth.v1.1-std.upf',
 'Cu.paw.z_19.atompaw.jth.v1.1-std.upf',
 'Zn.paw.z_12.atompaw.jth.v1.1-std.upf',
 'Ga.paw.z_13.atompaw.jth.v1.1-std.upf',
 'Ge.paw.z_14.atompaw.jth.v1.1-std.upf',
 'As.paw.z_15.atompaw.jth.v1.1-std.upf',
 'Se.paw.z_6.atompaw.jth.v1.1-std.upf', 
 'Br.paw.z_7.atompaw.jth.v1.1-std.upf', 
 'Kr.paw.z_8.atompaw.jth.v1.1-std.upf', 
 'Rb.paw.z_9.atompaw.jth.v1.1-std.upf', 
 'Sr.paw.z_10.atompaw.jth.v1.1-std.upf',
 'Y.paw.z_11.atompaw.jth.v1.1-std.upf', 
 'Zr.paw.z_12.atompaw.jth.v1.1-std.upf',
 'Nb.paw.z_13.atompaw.jth.v1.1-std.upf',
 'Mo.paw.z_14.atompaw.jth.v1.1-std.upf',
 'Tc.paw.z_15.atompaw.jth.v1.1-std.upf',
 'Ru.paw.z_16.atompaw.jth.v1.1-std.upf',
 'Rh.paw.z_17.atompaw.jth.v1.1-std.upf',
 'Pd.paw.z_18.atompaw.jth.v1.1-std.upf',
 'Ag.paw.z_11.atompaw.jth.v1.1-std.upf',
 'Cd.paw.z_12.atompaw.jth.v1.1-std.upf',
 'In.paw.z_13.atompaw.jth.v1.1-std.upf',
 'Sn.paw.z_14.atompaw.jth.v1.1-std.upf',
 'Sb.paw.z_15.atompaw.jth.v1.1-std.upf',
 'Te.paw.z_6.atompaw.jth.v1.1-std.upf',
 'I.paw.z_7.atompaw.jth.v1.1-std.upf',
 'Xe.paw.z_8.atompaw.jth.v1.1-std.upf',
 'Cs.paw.z_9.atompaw.jth.v1.1-std.upf',
 'Ba.paw.z_10.atompaw.jth.v1.1-std.upf',
 'Hf.paw.z_12.atompaw.jth.v1.1-std.upf',
 'Ta.paw.z_13.atompaw.jth.v1.1-std.upf',
 'W.paw.z_14.atompaw.jth.v1.1-std.upf', 
 'Re.paw.z_15.atompaw.jth.v1.1-std.upf',
 'Os.paw.z_8.atompaw.jth.v1.1-std.upf',
 'Ir.paw.z_15.atompaw.jth.v1.1-std.upf',  
 'Pt.paw.z_10.atompaw.jth.v1.1-std.upf', 
 'Au.paw.z_11.atompaw.jth.v1.1-std.upf', 
 'Hg.paw.z_12.atompaw.jth.v1.1-std.upf',
 'Tl.paw.z_13.atompaw.jth.v1.1-std.upf',
 'Pb.paw.z_14.atompaw.jth.v1.1-std.upf',
 'Bi.paw.z_15.atompaw.jth.v1.1-std.upf',
 'Po.paw.z_6.atompaw.jth.v1.1-std.upf', 
 'At.paw.z_7.atompaw.jth.v1.1-std.upf',
 'Rn.paw.z_8.atompaw.jth.v1.1-std.upf', 
 #'La.paw.z_11.atompaw.jth.v1.1-std.upf',
 #'Ce.paw.z_12.atompaw.jth.v1.1-std.upf',
 #'Pr.paw.z_13.atompaw.jth.v1.1-std.upf',
 #'Nd.paw.z_14.atompaw.jth.v1.1-std.upf',
 #'Pm.paw.z_15.atompaw.jth.v1.1-std.upf',
 #'Sm.paw.z_16.atompaw.jth.v1.1-std.upf',
 #'Eu.paw.z_17.atompaw.jth.v1.1-std.upf',
 #'Gd.paw.z_18.atompaw.jth.v1.1-std.upf',
 #'Tb.paw.z_19.atompaw.jth.v1.1-std.upf',
 #'Dy.paw.z_20.atompaw.jth.v1.1-std.upf',
 #'Ho.paw.z_21.atompaw.jth.v1.1-std.upf',
 #'Er.paw.z_22.atompaw.jth.v1.1-std.upf',
 #'Tm.paw.z_23.atompaw.jth.v1.1-std.upf',
 #'Yb.paw.z_24.atompaw.jth.v1.1-std.upf',
 #'Lu.paw.z_25.atompaw.jth.v1.1-std.upf',
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
