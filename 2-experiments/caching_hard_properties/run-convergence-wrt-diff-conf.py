# verdi run run-convergence-wrt-diff-conf.py

import os

pseudos = [
    "PAW-PSL1.0.0-high/Si.paw.z_4.ld1.psl.v1.0.0-high.upf",
    "PAW-PSL1.0.0-high/Cs.paw.z_9.ld1.psl.v1.0.0-high.upf",
    "PAW-PSL1.0.0-high/Cu.paw.z_19.ld1.psl.v1.0.0-high.upf",
    "PAW-PSL1.0.0-high/Fe.paw.z_16.ld1.psl.v1.0.0-high.upf",
    "PAW-PSL1.0.0-high/Te.paw.z_16.ld1.psl.v1.0.0-high.upf",
    "PAW-PSL1.0.0-high/Hf.paw.z_26.ld1.psl.v1.0.0-high.spfn.upf",
    "PAW-PSL1.0.0-high/Sn.paw.z_14.ld1.psl.v1.0.0-high.upf",
    "PAW-PSL1.0.0-high/Ar.paw.z_8.ld1.psl.v1.0.0-high.upf",
    "PAW-PSL1.0.0-high/U.paw.z_14.ld1.psl.v1.0.0-high.upf",
    "NC-DOJOv4-standard/Si.nc.z_4.oncvpsp3.dojo.v0.4.1-std.upf",
    "NC-DOJOv4-standard/Cs.nc.z_9.oncvpsp3.dojo.v0.4.1-std.upf",
    "NC-DOJOv4-standard/Cu.nc.z_19.oncvpsp3.dojo.v0.4.1-std.upf",
    "NC-DOJOv4-standard/Fe.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf",
    "NC-DOJOv4-standard/Te.nc.z_16.oncvpsp3.dojo.v0.4.1-std.upf",
    "NC-DOJOv4-standard/Hf.nc.z_12.oncvpsp3.dojo.v0.4.1-std.upf",
    "NC-DOJOv4-standard/Sn.nc.z_14.oncvpsp3.dojo.v0.4.1-std.upf",
    "NC-DOJOv4-standard/Ar.nc.z_8.oncvpsp3.dojo.v0.4.1-std.upf",
]

criteria = "efficiency"
# conf = "Diamond"
conf = None
base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe"

for pseudo in pseudos:
    pseudo_path = os.path.join(base_path, pseudo)
    if conf is not None:
        command = f"aiida-sssp-workflow launch --configuration {conf} --property convergence.phonon_frequencies --property convergence.bands --pw-code pw-7.0@eiger-mc-mr32-mem --ph-code ph-7.0@eiger-mc-mr32-mem --protocol acwf --cutoff-control standard --criteria {criteria} --withmpi True --num-mpiprocs 128 --npool 8  -- {pseudo_path}"
    else:
        command = f"aiida-sssp-workflow launch  --property convergence.phonon_frequencies --property convergence.bands --pw-code pw-7.0@eiger-mc-mr32-mem --ph-code ph-7.0@eiger-mc-mr32-mem --protocol acwf --cutoff-control standard --criteria {criteria} --withmpi True --num-mpiprocs 128 --npool 8  -- {pseudo_path}"
    os.system(command)
    # print(command)
    print(f"Launched {pseudo}")
