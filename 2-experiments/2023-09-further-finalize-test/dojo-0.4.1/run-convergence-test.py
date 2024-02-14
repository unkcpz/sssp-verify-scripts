# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "dojo-0.4.1-std"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.4-standard'


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
 #'H.nc.z_1.oncvpsp3.dojo.v0.4.1-std.upf',
 #'He.nc.z_2.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Li.nc.z_3.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Be.nc.z_4.oncvpsp3.dojo.v0.4.1-std.upf',
 #'B.nc.z_3.oncvpsp3.dojo.v0.4.1-std.upf',
 #'C.nc.z_4.oncvpsp3.dojo.v0.4.1-std.upf',
 #'N.nc.z_5.oncvpsp3.dojo.v0.4.1-std.upf',
 #'O.nc.z_6.oncvpsp3.dojo.v0.4.1-std.upf',
 #'F.nc.z_7.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ne.nc.z_8.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Na.nc.z_9.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Mg.nc.z_10.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Al.nc.z_3.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Si.nc.z_4.oncvpsp3.dojo.v0.4.1-std.upf',
 #'P.nc.z_5.oncvpsp3.dojo.v0.4.1-std.upf',
 'S.nc.z_6.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Cl.nc.z_7.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ar.nc.z_8.oncvpsp3.dojo.v0.4.1-std.upf',
 #'K.nc.z_9.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ca.nc.z_10.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Sc.nc.z_11.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ti.nc.z_12.oncvpsp3.dojo.v0.4.1-std.upf',
 #'V.nc.z_13.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Cr.nc.z_14.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Mn.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Fe.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Co.nc.z_17.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ni.nc.z_18.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Cu.nc.z_19.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Zn.nc.z_20.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ga.nc.z_13.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ge.nc.z_14.oncvpsp3.dojo.v0.4.1-std.upf',
 #'As.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Se.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Br.nc.z_7.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Kr.nc.z_8.oncvpsp3.dojo.v0.4.1-std.upf',
 'Rb.nc.z_9.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Sr.nc.z_10.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Y.nc.z_11.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Zr.nc.z_12.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Nb.nc.z_13.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Mo.nc.z_14.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Tc.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ru.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Rh.nc.z_17.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Pd.nc.z_18.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ag.nc.z_19.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Cd.nc.z_20.oncvpsp3.dojo.v0.4.1-std.upf',
 #'In.nc.z_13.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Sn.nc.z_14.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Sb.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf',
 'Te.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf',
 'I.nc.z_7.oncvpsp3.dojo.v0.4.1-std.upf',
 'Xe.nc.z_8.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Cs.nc.z_9.oncvpsp3.dojo.v0.4.1-std.upf',
 'Ba.nc.z_10.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Hf.nc.z_12.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ta.nc.z_13.oncvpsp3.dojo.v0.4.1-std.upf',
 #'W.nc.z_14.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Re.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Os.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Ir.nc.z_17.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Pt.nc.z_18.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Au.nc.z_19.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Hg.nc.z_20.oncvpsp3.dojo.v0.4.1-std.upf',
 'Tl.nc.z_13.oncvpsp3.dojo.v0.4.1-std.upf',
 'Pb.nc.z_14.oncvpsp3.dojo.v0.4.1-std.upf',
 'Bi.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf',
 'Po.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf',
 'Rn.nc.z_18.oncvpsp3.dojo.v0.4.1-std.upf',
 #'La.nc.z_11.oncvpsp3.dojo.v0.4.1-std.upf',
 #'Lu.nc.z_25.oncvpsp3.dojo.v0.4.1-std.upf',
    ]


    computer = 'eiger-mc-mr32-mem'
    #computer = 'daint-mc-mrcloud-mem'
    mpiprocs = 128 # 128 for eiger
    npool = 16 # 16 for eiger
    base_path = lib_path

    for pseudo in pseudos:
        pseudo_path = os.path.join(base_path, pseudo)
        command = f"aiida-sssp-workflow launch --property convergence --pw-code pw-7.0@{computer} --ph-code ph-7.0@{computer} --protocol acwf --cutoff-control standard --criteria {criteria} --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment}  -- {pseudo_path}"
        os.system(command)
        #print(command)
        print(f"Launched {pseudo}")
