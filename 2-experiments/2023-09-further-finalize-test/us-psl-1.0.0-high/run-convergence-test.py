# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "US-PSL-1.0.0-high"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-PSL1.0.0-high'


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
 'H.us.z_1.ld1.psl.v1.0.0-high.upf',  
 'He.us.z_2.ld1.psl.v1.0.0-high.upf', 
 'Li.us.z_3.ld1.psl.v1.0.0-high.upf', 
 'Be.us.z_4.ld1.psl.v1.0.0-high.upf', 
 'B.us.z_3.ld1.psl.v1.0.0-high.upf',  
 'C.us.z_4.ld1.psl.v1.0.0-high.upf',  
 'N.us.z_5.ld1.psl.v1.0.0-high.upf',  
 'O.us.z_6.ld1.psl.v1.0.0-high.upf',  
 'F.us.z_7.ld1.psl.v1.0.0-high.upf',  
 'Ne.us.z_8.ld1.psl.v1.0.0-high.upf', 
 'Na.us.z_9.ld1.psl.v1.0.0-high.upf', 
 'Mg.us.z_10.ld1.psl.v1.0.0-high.upf',
 'Al.us.z_3.ld1.psl.v1.0.0-high.upf', 
 'Si.us.z_4.ld1.psl.v1.0.0-high.upf', 
 'P.us.z_5.ld1.psl.v1.0.0-high.upf',  
 'S.us.z_6.ld1.psl.v1.0.0-high.upf',  
 'Cl.us.z_7.ld1.psl.v1.0.0-high.upf',      
 'Ar.us.z_8.ld1.psl.v1.0.0-high.upf',       
 'K.us.z_9.ld1.psl.v1.0.0-high.upf',  
 'Ca.us.z_10.ld1.psl.v1.0.0-high.upf',
 'Sc.us.z_11.ld1.psl.v1.0.0-high.upf',
 'Ti.us.z_12.ld1.psl.v1.0.0-high.upf',
 'V.us.z_13.ld1.psl.v1.0.0-high.upf', 
 'Cr.us.z_14.ld1.psl.v1.0.0-high.upf',
 'Mn.us.z_15.ld1.psl.v1.0.0-high.upf',
 'Fe.us.z_16.ld1.psl.v1.0.0-high.upf',
 'Co.us.z_17.ld1.psl.v1.0.0-high.upf',
 'Ni.us.z_18.ld1.psl.v1.0.0-high.upf',      
 'Cu.us.z_19.ld1.psl.v1.0.0-high.upf',     
 'Zn.us.z_20.ld1.psl.v1.0.0-high.upf',
 'Ga.us.z_13.ld1.psl.v1.0.0-high.upf',
 'Ge.us.z_14.ld1.psl.v1.0.0-high.upf',
 'As.us.z_15.ld1.psl.v1.0.0-high.upf',
 'Se.us.z_16.ld1.psl.v1.0.0-high.upf',
 'Br.us.z_17.ld1.psl.v1.0.0-high.upf',
 'Kr.us.z_18.ld1.psl.v1.0.0-high.upf',
 'Rb.us.z_9.ld1.psl.v1.0.0-high.upf', 
 'Sr.us.z_10.ld1.psl.v1.0.0-high.upf',
 'Y.us.z_11.ld1.psl.v1.0.0-high.upf',           
 'Zr.us.z_12.ld1.psl.v1.0.0-high.upf',
 'Nb.us.z_13.ld1.psl.v1.0.0-high.upf',
 'Mo.us.z_14.ld1.psl.v1.0.0-high.upf',
 'Tc.us.z_15.ld1.psl.v1.0.0-high.upf',
 'Ru.us.z_16.ld1.psl.v1.0.0-high.upf',
 'Rh.us.z_17.ld1.psl.v1.0.0-high.upf',
 'Pd.us.z_18.ld1.psl.v1.0.0-high.upf',
 'Ag.us.z_19.ld1.psl.v1.0.0-high.upf',
 'Cd.us.z_20.ld1.psl.v1.0.0-high.upf',
 'In.us.z_13.ld1.psl.v1.0.0-high.upf',
 'Sn.us.z_14.ld1.psl.v1.0.0-high.upf',
 'Sb.us.z_15.ld1.psl.v1.0.0-high.upf',
 'Te.us.z_16.ld1.psl.v1.0.0-high.upf',
 'I.us.z_17.ld1.psl.v1.0.0-high.upf',
 'Xe.us.z_18.ld1.psl.v1.0.0-high.upf',
 'Cs.us.z_9.ld1.psl.v1.0.0-high.upf',
 'Ba.us.z_10.ld1.psl.v1.0.0-high.upf',
 'Hf.us.z_36.ld1.psl.v1.0.0-high.spdfn.upf',
 'Hf.us.z_26.ld1.psl.v1.0.0-high.spfn.upf',
 'Ta.us.z_27.ld1.psl.v1.0.0-high.upf',
 'W.us.z_28.ld1.psl.v1.0.0-high.upf',
 'Re.us.z_29.ld1.psl.v1.0.0-high.upf',
 'Os.us.z_30.ld1.psl.v1.0.0-high.upf',
 'Ir.us.z_31.ld1.psl.v1.0.0-high.upf',
 'Pt.us.z_32.ld1.psl.v1.0.0-high.upf',
 'Au.us.z_33.ld1.psl.v1.0.0-high.upf',
 'Hg.us.z_20.ld1.psl.v1.0.0-high.upf',
 'Tl.us.z_13.ld1.psl.v1.0.0-high.upf',
 'Pb.us.z_14.ld1.psl.v1.0.0-high.upf',
 'Bi.us.z_15.ld1.psl.v1.0.0-high.upf',
 'Po.us.z_16.ld1.psl.v1.0.0-high.upf',
 'At.us.z_17.ld1.psl.v1.0.0-high.upf',
 'Rn.us.z_18.ld1.psl.v1.0.0-high.upf',
 'Fr.us.z_19.ld1.psl.v1.0.0-high.upf',
 'Ra.us.z_20.ld1.psl.v1.0.0-high.upf',
 #'La.us.z_11.ld1.psl.v1.0.0-high.spfn.upf',
 #'La.us.z_21.ld1.psl.v1.0.0-high.spdfn.upf',
 #'Ce.us.z_22.ld1.psl.v1.0.0-high.upf',
 #'Pr.us.z_23.ld1.psl.v1.0.0-high.upf',
 #'Nd.us.z_24.ld1.psl.v1.0.0-high.upf',
 #'Pm.us.z_25.ld1.psl.v1.0.0-high.upf',
 #'Sm.us.z_26.ld1.psl.v1.0.0-high.upf',
 #'Eu.us.z_27.ld1.psl.v1.0.0-high.upf',
 #'Gd.us.z_28.ld1.psl.v1.0.0-high.upf',
 #'Tb.us.z_29.ld1.psl.v1.0.0-high.upf',
 #'Dy.us.z_30.ld1.psl.v1.0.0-high.upf',
 #'Ho.us.z_31.ld1.psl.v1.0.0-high.upf',
 #'Er.us.z_32.ld1.psl.v1.0.0-high.upf',
 #'Tm.us.z_33.ld1.psl.v1.0.0-high.upf',
 #'Yb.us.z_34.ld1.psl.v1.0.0-high.upf',
 #'Lu.us.z_35.ld1.psl.v1.0.0-high.upf',
 #'Th.us.z_12.ld1.psl.v1.0.0-high.upf',
 #'Pa.us.z_13.ld1.psl.v1.0.0-high.upf',
 #'U.us.z_14.ld1.psl.v1.0.0-high.upf',
 #'Np.us.z_15.ld1.psl.v1.0.0-high.upf',
 #'Pu.us.z_16.ld1.psl.v1.0.0-high.upf',
 #'Am.us.z_17.ld1.psl.v1.0.0-high.upf',
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
