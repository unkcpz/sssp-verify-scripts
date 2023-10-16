# verdi run run-convergence-wrt-diff-conf.py

import os

pseudos = [
    #"PAW-PSL1.0.0-high/Fe.paw.z_16.ld1.psl.v1.0.0-high.upf",
    #"PAW-PSL1.0.0-high/Cr.paw.z_14.ld1.psl.v1.0.0-high.upf",
    "PAW-PSL1.0.0-high/O.paw.z_6.ld1.psl.v1.0.0-high.upf",
    "US-PSL1.0.0-low/O.us.z_6.ld1.psl.v1.0.0-low.upf",
    #"NC-DOJOv4-standard/Fe.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf",
    #"NC-DOJOv4-standard/Cr.nc.z_14.oncvpsp3.dojo.v0.4.1-std.upf",
    "NC-DOJOv4-standard/O.nc.z_6.oncvpsp3.dojo.v0.4.1-std.upf",
]

#criteria = "efficiency"
criteria = "precision"
conf = 'GS'
comment = "mag-on"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud-mem'
mpiprocs = 32 # 128 for eiger
npool = 4 # 8 for eiger
base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe"

for pseudo in pseudos:
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --configuration {conf} --property convergence.cohesive_energy --pw-code pw-7.0@{computer} --ph-code ph-7.0@{computer} --protocol acwf --cutoff-control standard --criteria {criteria} --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment}  -- {pseudo_path}"
    os.system(command)
    # print(command)
    print(f"Launched {pseudo}")
