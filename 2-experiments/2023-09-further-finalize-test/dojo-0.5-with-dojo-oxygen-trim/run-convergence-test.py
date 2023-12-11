# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "dojo-0.5-standard-trim"

lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.5-standard-trim'


if mode == "pseudos_generate":
    pd_elements_lst = [
        "Bi", "Cd", "Cu", "Pb", "Pd", "Po", "Rn", "Te", "Tl", "Xe", "Zn",
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
    pseudos = ['Bi.nc.z_15.oncvpsp4.dojo.v0.5.0-std-trim.upf',
 'Cd.nc.z_20.oncvpsp3.dojo.v0.5.0-std-trim.upf',
 'Cu.nc.z_19.oncvpsp3.dojo.v0.5.0-std-trim.upf',
 'Pb.nc.z_14.oncvpsp4.dojo.v0.5.0-std-trim.upf',
 'Pd.nc.z_18.oncvpsp3.dojo.v0.5.0-std-trim.upf',
 'Po.nc.z_16.oncvpsp4.dojo.v0.5.0-std-trim.upf',
 'Rn.nc.z_18.oncvpsp4.dojo.v0.5.0-std-trim.upf',
 'Te.nc.z_16.oncvpsp4.dojo.v0.5.0-std-trim.upf',
 'Tl.nc.z_13.oncvpsp4.dojo.v0.5.0-std-trim.upf',
 'Xe.nc.z_8.oncvpsp4.dojo.v0.5.0-std-trim.upf',
 'Zn.nc.z_20.oncvpsp3.dojo.v0.5.0-std-trim.upf',
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
