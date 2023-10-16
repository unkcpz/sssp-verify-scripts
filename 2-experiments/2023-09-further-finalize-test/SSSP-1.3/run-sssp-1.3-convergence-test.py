# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "SSSP-1.3-precision"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/MIX-SSSP-precision-1.3.0-recollected'


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

pseudos = [
    'H.nc.z_1.oncvpsp3.sg15.v1.0.upf',
    'He.nc.z_2.oncvpsp3.sg15.v1.0.upf',
    'Li.us.z_3.uspp.gbrv.v1.4.upf', 
    'Be.nc.z_4.oncvpsp3.sg15.v1.0.upf',
    'B.us.z_3.uspp.gbrv.v1.4.upf',       
    'C.paw.z_4.ld1.psl.v1.0.0-high.upf', 
    'N.nc.z_5.oncvpsp3.dojo.v0.4.1-std.upf',
    'O.paw.z_6.ld1.psl.v0.1.upf',   
    'F.nc.z_7.oncvpsp3.dojo.v0.4.1-std.upf',
    'Ne.paw.z_8.ld1.psl.v1.0.0-high.upf',
    'Na.paw.z_9.ld1.psl.v1.0.0-low.upf',  
    'Mg.us.z_10.uspp.gbrv.v1.4.upf',     
    'Al.paw.z_3.ld1.psl.v1.0.0-high.upf',     
    'Si.us.z_4.ld1.psl.v1.0.0-high.upf',        
    'P.us.z_5.ld1.psl.v1.0.0-high.upf',         
    'S.us.z_6.uspp.gbrv.v1.4.upf',              
    'Cl.us.z_7.ld1.psl.v1.0.0-high.upf',        
    'Ar.paw.z_8.ld1.psl.v1.0.0-high.upf',       
    'K.paw.z_9.ld1.psl.v1.0.0-high.upf',        
    'Ca.us.z_10.uspp.gbrv.v1.upf',              
    'Sc.paw.z_11.ld1.psl.v0.2.3.upf',           
    'Ti.us.z_12.uspp.gbrv.v1.4.upf',            
    'V.us.z_13.uspp.gbrv.v1.4.upf',             
    'Cr.us.z_14.uspp.gbrv.v1.5.upf',            
    'Mn.us.z_15.uspp.gbrv.v1.5.upf',            
    'Fe.paw.z_16.ld1.psl.v0.2.1.upf',           
    'Co.us.z_17.uspp.gbrv.v1.2.upf',            
    'Ni.us.z_18.uspp.gbrv.v1.4.upf',            
    'Cu.paw.z_11.ld1.psl.v1.0.0-low.upf',
    'Zn.us.z_20.uspp.gbrv.v1.upf',       
    'Ga.paw.z_13.ld1.psl.v1.0.0-high.upf',
    'Ge.us.z_14.uspp.gbrv.v1.4.upf',     
    'As.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf',
    'Se.us.z_6.uspp.gbrv.v1.upf',         
    'Br.us.z_7.uspp.gbrv.v1.4.upf',          
    'Kr.paw.z_18.ld1.psl.v1.0.0-high.upf',
    'Rb.nc.z_9.oncvpsp3.sg15.v1.0.upf',       
    'Sr.us.z_10.uspp.gbrv.v1.upf',       
    'Y.us.z_11.uspp.gbrv.v1.4.upf',      
    'Zr.us.z_12.uspp.gbrv.v1.upf',       
    'Nb.paw.z_13.ld1.psl.v0.3.0.upf',    
    'Mo.nc.z_14.oncvpsp3.sg15.v1.0.upf', 
    'Tc.nc.z_15.oncvpsp3.sg15.v1.0.upf', 
    'Ru.nc.z_16.oncvpsp3.sg15.v1.0.upf',         
    'Rh.nc.z_17.oncvpsp3.sg15.v1.0.upf',
    'Pd.nc.z_18.oncvpsp3.sg15.v1.0.upf',
    'Ag.nc.z_19.oncvpsp3.sg15.v1.0.upf',
    'Cd.paw.z_20.ld1.psl.v1.0.0-high.upf',
    'In.us.z_13.ld1.psl.v0.2.2.upf',
    'Sn.us.z_14.uspp.gbrv.v1.4.upf',
    'Sb.us.z_15.uspp.gbrv.v1.4.upf',
    'Te.us.z_6.ld1.psl.v1.0.0-low.upf',
    'I.nc.z_17.oncvpsp4.sg15.v0.upf',
    'Xe.paw.z_18.ld1.psl.v1.0.0-high.upf',
    'Cs.nc.z_9.oncvpsp3.dojo.v0.4.1-str.upf',
    'Ba.nc.z_10.oncvpsp4.dojo.v0.5.0.upf',
    'Hf.nc.z_12.oncvpsp3.dojo.v0.4.1-std.upf',
    'Ta.us.z_13.uspp.gbrv.v1.upf',
    'W.us.z_14.uspp.gbrv.v1.2.upf',
    'Re.us.z_15.uspp.gbrv.v1.2.upf',
    'Os.us.z_16.uspp.gbrv.v1.2.upf',
    'Ir.us.z_31.ld1.psl.v1.0.0-high.upf',
    'Pt.us.z_32.ld1.psl.v1.0.0-high.upf',
    'Au.nc.z_19.oncvpsp3.sg15.v1.0.upf',
    'Hg.us.z_12.uspp.gbrv.v1.upf',
    'Tl.us.z_13.uspp.gbrv.v1.2.upf',
    'Pb.paw.z_14.ld1.psl.v0.2.2.upf',
    'Bi.us.z_15.uspp.gbrv.v1.upf',
    'Po.us.z_16.ld1.psl.v1.0.0-high.upf',
    'Rn.paw.z_18.ld1.psl.v1.0.0-high.upf',
    #'La.paw.z_11.atompaw.wentzcovitch.v1.2.upf',
    #'Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf',
    #'Pr.paw.z_13.atompaw.wentzcovitch.v1.2.upf',
    #'Nd.paw.z_14.atompaw.wentzcovitch.v1.2.upf',
    #'Pm.paw.z_15.atompaw.wentzcovitch.v1.2.upf',
    #'Sm.paw.z_16.atompaw.wentzcovitch.v1.2.upf',
    #'Eu.paw.z_17.atompaw.wentzcovitch.v1.2.upf',
    #'Gd.paw.z_18.atompaw.wentzcovitch.v1.2.upf',
    #'Tb.paw.z_19.atompaw.wentzcovitch.v1.2.upf',
    #'Dy.paw.z_20.atompaw.wentzcovitch.v1.2.upf',
    #'Ho.paw.z_21.atompaw.wentzcovitch.v1.2.upf',
    #'Er.paw.z_22.atompaw.wentzcovitch.v1.2.upf',
    #'Tm.paw.z_23.atompaw.wentzcovitch.v1.2.upf',
    #'Yb.paw.z_24.atompaw.wentzcovitch.v1.2.upf',
    #'Lu.paw.z_25.atompaw.wentzcovitch.v1.2.upf',
    #'Th.paw.z_12.ld1.uni-marburg.v0.upf',
    #'Pa.paw.z_13.ld1.uni-marburg.v0.upf',
    #'U.paw.z_14.ld1.uni-marburg.v0.upf',
    #'Np.paw.z_15.ld1.uni-marburg.v0.upf',
    #'Pu.paw.z_16.ld1.uni-marburg.v0.upf',
    #'Am.paw.z_17.ld1.uni-marburg.v0.upf',
    #'Cm.paw.z_18.ld1.uni-marburg.v0.upf',
    #'Bk.paw.z_19.ld1.uni-marburg.v0.upf',
    #'Cf.paw.z_20.ld1.uni-marburg.v0.upf',
    #'Es.paw.z_21.ld1.uni-marburg.v0.upf',
    #'Fm.paw.z_22.ld1.uni-marburg.v0.upf',
    #'Md.paw.z_23.ld1.uni-marburg.v0.upf',
    #'No.paw.z_24.ld1.uni-marburg.v0.upf',
    #'Lr.paw.z_25.ld1.uni-marburg.v0.upf',
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
