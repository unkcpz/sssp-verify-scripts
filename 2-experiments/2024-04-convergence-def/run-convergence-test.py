# verdi run run-sssp-1.3-convergence-test.py
import os
from pathlib import Path
import pprint

mode = "run" # or "pseudos_generate"
#mode = "pseudos_generate"
#criteria = "efficiency"
criteria = "precision"
comment = "convergence-def"

gbrv_lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-GBRV-1.x'
jth_lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-JTH1.1-standard'
dojo_lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.4-standard'
psl_paw_lib_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-PSL1.0.0-high'


# Hg, Ga, N, Cs, Mn
pseudos = [
    # GBRV-1.x
    f"{gbrv_lib_path}/Hg.us.z_12.uspp.gbrv.v1.upf",
    f"{gbrv_lib_path}/Ga.us.z_19.uspp.gbrv.v1.4.upf",
    f"{gbrv_lib_path}/N.us.z_5.uspp.gbrv.v1.2.upf",
    f"{gbrv_lib_path}/Cs.us.z_9.uspp.gbrv.v1.upf",
    f"{gbrv_lib_path}/Mn.us.z_15.uspp.gbrv.v1.5.upf",
    # JTH
    f"{jth_lib_path}/Hg.paw.z_12.atompaw.jth.v1.1-std.upf",
    f"{jth_lib_path}/Ga.paw.z_13.atompaw.jth.v1.1-std.upf",
    f"{jth_lib_path}/N.paw.z_5.atompaw.jth.v1.1-std.upf",
    f"{jth_lib_path}/Cs.paw.z_9.atompaw.jth.v1.1-std.upf",
    f"{jth_lib_path}/Mn.paw.z_15.atompaw.jth.v1.1-std.upf",
    # DOJO
    f"{dojo_lib_path}/Hg.nc.z_20.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/Ga.nc.z_13.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/N.nc.z_5.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/Cs.nc.z_9.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/Mn.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf",
    # PSL-PAW
    f"{psl_paw_lib_path}/Hg.paw.z_20.ld1.psl.v1.0.0-high.upf",
    f"{psl_paw_lib_path}/Ga.paw.z_13.ld1.psl.v1.0.0-high.upf",
    f"{psl_paw_lib_path}/N.paw.z_5.ld1.psl.v1.0.0-high.upf",
    f"{psl_paw_lib_path}/Cs.paw.z_9.ld1.psl.v1.0.0-high.upf",
    f"{psl_paw_lib_path}/Mn.paw.z_15.ld1.psl.v1.0.0-high.upf",
]

computer = 'daint-mr32-mem'
mpiprocs = 36
npool = 4
conf = 'BCC'

for pseudo_path in pseudos:
    pseudo = Path(pseudo_path).name
    command = f"aiida-sssp-workflow launch --property convergence --configuration {conf} --pw-code pw-7.0@{computer} --ph-code ph-7.0@{computer} --protocol acwf --cutoff-control refined --criteria {criteria} --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment}  -- {pseudo_path}"
    os.system(command)
    #print(command)
    # check the path exists
    print(f"Launched {pseudo}")
