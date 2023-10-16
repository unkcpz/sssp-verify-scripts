# verdi run run-convergence-wrt-diff-conf.py

import os

pseudos = [
"N.nc.z_5.oncvpsp3.dojo.v0.4.1-std.upf",
"N.nc.z_5.oncvpsp4.sg15.v0.upf",
"N.nc.z_5.oncvpsp4.spms.v1.upf",
"N.paw.z_5.atompaw.jth.v1.1-std.upf",
"N.paw.z_5.atompaw.jth.v1.1-str.upf",
"N.paw.z_5.ld1.psl.v0.1.upf",
"N.paw.z_5.ld1.psl.v1.0.0-high.upf",
"N.us.z_5.ld1.psl.v0.1.upf",
"N.us.z_5.ld1.psl.v1.0.0-high.upf",
"N.us.z_5.ld1.theose.v0.upf",
"N.us.z_5.uspp.gbrv.v1.2.upf",
]

#criteria = "efficiency"
criteria = "precision"
comment = "Nitrigen"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud-mem'
mpiprocs = 32 # 128 for eiger
npool = 4 # 8 for eiger
base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/_sssp_pbe/N/"

for pseudo in pseudos:
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property convergence --pw-code pw-7.0@{computer} --ph-code ph-7.0@{computer} --protocol acwf --cutoff-control standard --criteria {criteria} --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment}  -- {pseudo_path}"
    os.system(command)
    # print(command)
    print(f"Launched {pseudo}")
