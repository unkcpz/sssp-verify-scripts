# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "PAW-PSL-1.0.0-high"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-PSL1.0.0-high'


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
 'H.paw.z_1.ld1.psl.v1.0.0-high.upf',  
 'He.paw.z_2.ld1.psl.v1.0.0-high.upf', 
 'Li.paw.z_3.ld1.psl.v1.0.0-high.upf', 
 'Be.paw.z_4.ld1.psl.v1.0.0-high.upf', 
 'B.paw.z_3.ld1.psl.v1.0.0-high.upf',  
 'C.paw.z_4.ld1.psl.v1.0.0-high.upf',  
 'N.paw.z_5.ld1.psl.v1.0.0-high.upf',  
 'O.paw.z_6.ld1.psl.v1.0.0-high.upf',  
 'F.paw.z_7.ld1.psl.v1.0.0-high.upf',  
 'Ne.paw.z_8.ld1.psl.v1.0.0-high.upf', 
 'Na.paw.z_9.ld1.psl.v1.0.0-high.upf', 
 'Mg.paw.z_10.ld1.psl.v1.0.0-high.upf',
 'Al.paw.z_3.ld1.psl.v1.0.0-high.upf', 
 'Si.paw.z_4.ld1.psl.v1.0.0-high.upf', 
 'P.paw.z_5.ld1.psl.v1.0.0-high.upf',  
 'S.paw.z_6.ld1.psl.v1.0.0-high.upf',  
 'Cl.paw.z_7.ld1.psl.v1.0.0-high.upf',       
 'Ar.paw.z_8.ld1.psl.v1.0.0-high.upf',      
 'K.paw.z_9.ld1.psl.v1.0.0-high.upf',  
 'Ca.paw.z_10.ld1.psl.v1.0.0-high.upf',
 'Sc.paw.z_11.ld1.psl.v1.0.0-high.upf',
 'Ti.paw.z_12.ld1.psl.v1.0.0-high.upf',
 'V.paw.z_13.ld1.psl.v1.0.0-high.upf', 
 'Cr.paw.z_14.ld1.psl.v1.0.0-high.upf',
 #'Mn.paw.z_15.ld1.psl.v1.0.0-high.upf',
 #'Fe.paw.z_16.ld1.psl.v1.0.0-high.upf',
 #'Co.paw.z_17.ld1.psl.v1.0.0-high.upf',
 'Ni.paw.z_18.ld1.psl.v1.0.0-high.upf',     
 'Cu.paw.z_19.ld1.psl.v1.0.0-high.upf',      
 'Zn.paw.z_20.ld1.psl.v1.0.0-high.upf',
 'Ga.paw.z_13.ld1.psl.v1.0.0-high.upf',
 'Ge.paw.z_14.ld1.psl.v1.0.0-high.upf',
 'As.paw.z_15.ld1.psl.v1.0.0-high.upf',
 'Se.paw.z_16.ld1.psl.v1.0.0-high.upf',
 'Br.paw.z_17.ld1.psl.v1.0.0-high.upf',
 'Kr.paw.z_18.ld1.psl.v1.0.0-high.upf',
 'Rb.paw.z_9.ld1.psl.v1.0.0-high.upf', 
 'Sr.paw.z_10.ld1.psl.v1.0.0-high.upf',
 'Y.paw.z_11.ld1.psl.v1.0.0-high.upf',
 'Zr.paw.z_12.ld1.psl.v1.0.0-high.upf',
 'Nb.paw.z_13.ld1.psl.v1.0.0-high.upf',
 'Mo.paw.z_14.ld1.psl.v1.0.0-high.upf',
 'Tc.paw.z_15.ld1.psl.v1.0.0-high.upf',
 'Ru.paw.z_16.ld1.psl.v1.0.0-high.upf',
 'Rh.paw.z_17.ld1.psl.v1.0.0-high.upf',
 'Pd.paw.z_18.ld1.psl.v1.0.0-high.upf',
 'Ag.paw.z_19.ld1.psl.v1.0.0-high.upf',
 'Cd.paw.z_20.ld1.psl.v1.0.0-high.upf',
 'In.paw.z_13.ld1.psl.v1.0.0-high.upf',
 'Sn.paw.z_14.ld1.psl.v1.0.0-high.upf',
 'Sb.paw.z_15.ld1.psl.v1.0.0-high.upf',
 'Te.paw.z_16.ld1.psl.v1.0.0-high.upf',
 'I.paw.z_17.ld1.psl.v1.0.0-high.upf',
 'Xe.paw.z_18.ld1.psl.v1.0.0-high.upf',
 'Cs.paw.z_9.ld1.psl.v1.0.0-high.upf',
 'Ba.paw.z_10.ld1.psl.v1.0.0-high.upf',
 'Hf.paw.z_26.ld1.psl.v1.0.0-high.spfn.upf',
 'Hf.paw.z_36.ld1.psl.v1.0.0-high.spdfn.upf',
 'Ta.paw.z_27.ld1.psl.v1.0.0-high.upf',
 'W.paw.z_28.ld1.psl.v1.0.0-high.upf',
 'Re.paw.z_29.ld1.psl.v1.0.0-high.upf',
 'Os.paw.z_30.ld1.psl.v1.0.0-high.upf',
 'Ir.paw.z_31.ld1.psl.v1.0.0-high.upf',
 'Pt.paw.z_32.ld1.psl.v1.0.0-high.upf',
 'Au.paw.z_33.ld1.psl.v1.0.0-high.upf',
 'Hg.paw.z_20.ld1.psl.v1.0.0-high.upf',
 'Tl.paw.z_13.ld1.psl.v1.0.0-high.upf',
 'Pb.paw.z_14.ld1.psl.v1.0.0-high.upf',
 'Bi.paw.z_15.ld1.psl.v1.0.0-high.upf',
 'Po.paw.z_16.ld1.psl.v1.0.0-high.upf',
 'At.paw.z_17.ld1.psl.v1.0.0-high.upf',
 'Rn.paw.z_18.ld1.psl.v1.0.0-high.upf',
 'Fr.paw.z_19.ld1.psl.v1.0.0-high.upf',
 'Ra.paw.z_20.ld1.psl.v1.0.0-high.upf',
 #'La.paw.z_21.ld1.psl.v1.0.0-high.spdfn.upf',
 #'La.paw.z_11.ld1.psl.v1.0.0-high.spfn.upf',
 #'Ce.paw.z_22.ld1.psl.v1.0.0-high.upf',
 #'Pr.paw.z_23.ld1.psl.v1.0.0-high.upf',
 #'Nd.paw.z_24.ld1.psl.v1.0.0-high.upf',
 #'Pm.paw.z_25.ld1.psl.v1.0.0-high.upf',
 #'Sm.paw.z_26.ld1.psl.v1.0.0-high.upf',
 #'Eu.paw.z_27.ld1.psl.v1.0.0-high.upf',
 #'Gd.paw.z_28.ld1.psl.v1.0.0-high.upf',
 #'Tb.paw.z_29.ld1.psl.v1.0.0-high.upf',
 #'Dy.paw.z_30.ld1.psl.v1.0.0-high.upf',
 #'Ho.paw.z_31.ld1.psl.v1.0.0-high.upf',
 #'Er.paw.z_32.ld1.psl.v1.0.0-high.upf',
 #'Tm.paw.z_33.ld1.psl.v1.0.0-high.upf',
 #'Yb.paw.z_34.ld1.psl.v1.0.0-high.upf',
 #'Lu.paw.z_35.ld1.psl.v1.0.0-high.upf',
 #'Th.paw.z_12.ld1.psl.v1.0.0-high.upf',
 #'Pa.paw.z_13.ld1.psl.v1.0.0-high.upf',
 #'U.paw.z_14.ld1.psl.v1.0.0-high.upf',
 #'Np.paw.z_15.ld1.psl.v1.0.0-high.upf',
 #'Pu.paw.z_16.ld1.psl.v1.0.0-high.upf',
 #'Am.paw.z_17.ld1.psl.v1.0.0-high.upf',
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
