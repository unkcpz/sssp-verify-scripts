# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "sg15"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-SG15-ONCVPSP4'


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
 'H.nc.z_1.oncvpsp4.sg15.v0.upf',
 'He.nc.z_2.oncvpsp4.sg15.v0.upf',
 'Li.nc.z_3.oncvpsp4.sg15.v0.upf',
 'Be.nc.z_4.oncvpsp4.sg15.v0.upf',
 'B.nc.z_3.oncvpsp4.sg15.v0.upf',
 'C.nc.z_4.oncvpsp4.sg15.v0.upf',
 'N.nc.z_5.oncvpsp4.sg15.v0.upf',
 'O.nc.z_6.oncvpsp4.sg15.v0.upf',
 'F.nc.z_7.oncvpsp4.sg15.v0.upf',
 'Ne.nc.z_8.oncvpsp4.sg15.v0.upf',
 'Na.nc.z_9.oncvpsp4.sg15.v0.upf',
 'Mg.nc.z_10.oncvpsp4.sg15.v0.upf',
 'Al.nc.z_11.oncvpsp4.sg15.v0.upf',
 'Si.nc.z_4.oncvpsp4.sg15.v0.upf',
 'P.nc.z_5.oncvpsp4.sg15.v0.upf',
 'S.nc.z_6.oncvpsp4.sg15.v0.upf',
 'Cl.nc.z_7.oncvpsp4.sg15.v0.upf',
 'K.nc.z_9.oncvpsp4.sg15.v0.upf',
 'Ca.nc.z_10.oncvpsp4.sg15.v0.upf',
 'Sc.nc.z_11.oncvpsp4.sg15.v0.upf',
 'Ti.nc.z_12.oncvpsp4.sg15.v0.upf',
 'V.nc.z_13.oncvpsp4.sg15.v0.upf',
 'Cr.nc.z_14.oncvpsp4.sg15.v0.upf',
 'Mn.nc.z_15.oncvpsp4.sg15.v0.upf',
 'Fe.nc.z_16.oncvpsp4.sg15.v0.upf',
 'Co.nc.z_17.oncvpsp4.sg15.v0.upf',
 'Ni.nc.z_18.oncvpsp4.sg15.v0.upf',
 'Cu.nc.z_19.oncvpsp4.sg15.v0.upf',
 'Zn.nc.z_20.oncvpsp4.sg15.v0.upf',
 'Ga.nc.z_13.oncvpsp4.sg15.v0.upf',
 'Ge.nc.z_14.oncvpsp4.sg15.v0.upf',
 'As.nc.z_5.oncvpsp4.sg15.v0.upf',
 'Se.nc.z_6.oncvpsp4.sg15.v0.upf',
 'Br.nc.z_7.oncvpsp4.sg15.v0.upf',
 'Kr.nc.z_8.oncvpsp4.sg15.v0.upf',
 'Rb.nc.z_9.oncvpsp4.sg15.v0.upf',
 'Sr.nc.z_10.oncvpsp4.sg15.v0.upf',
 'Y.nc.z_11.oncvpsp4.sg15.v0.upf',
 'Zr.nc.z_12.oncvpsp4.sg15.v0.upf',
 'Nb.nc.z_13.oncvpsp4.sg15.v0.upf',
 'Mo.nc.z_14.oncvpsp4.sg15.v0.upf',
 'Tc.nc.z_15.oncvpsp4.sg15.v0.upf',
 'Ru.nc.z_16.oncvpsp4.sg15.v0.upf',
 'Rh.nc.z_17.oncvpsp4.sg15.v0.upf',
 'Pd.nc.z_18.oncvpsp4.sg15.v0.upf',
 'Ag.nc.z_19.oncvpsp4.sg15.v0.upf',
 'Cd.nc.z_20.oncvpsp4.sg15.v0.upf',
 'In.nc.z_13.oncvpsp4.sg15.v0.upf',
 'Sn.nc.z_14.oncvpsp4.sg15.v0.upf',
 'Sb.nc.z_15.oncvpsp4.sg15.v0.upf',
 'Te.nc.z_16.oncvpsp4.sg15.v0.upf',
 'I.nc.z_17.oncvpsp4.sg15.v0.upf',
 'Xe.nc.z_18.oncvpsp4.sg15.v0.upf',
 'Cs.nc.z_9.oncvpsp4.sg15.v0.upf',
 'Ba.nc.z_10.oncvpsp4.sg15.v0.upf',
 'Hf.nc.z_26.oncvpsp4.sg15.v0.upf',
 'Ta.nc.z_27.oncvpsp4.sg15.v0.upf',
 'W.nc.z_28.oncvpsp4.sg15.v0.upf',
 'Re.nc.z_15.oncvpsp4.sg15.v0.upf',
 'Os.nc.z_16.oncvpsp4.sg15.v0.upf',
 'Ir.nc.z_17.oncvpsp4.sg15.v0.upf',
 'Pt.nc.z_18.oncvpsp4.sg15.v0.upf',
 'Au.nc.z_19.oncvpsp4.sg15.v0.upf',
 'Hg.nc.z_20.oncvpsp4.sg15.v0.upf',
 'Tl.nc.z_13.oncvpsp4.sg15.v0.upf',
 'Pb.nc.z_14.oncvpsp4.sg15.v0.upf',
 'Bi.nc.z_15.oncvpsp4.sg15.v0.upf',
 'La.nc.z_11.oncvpsp4.sg15.v0.upf',
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
