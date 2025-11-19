import os

pseudos = [
"O.nc.z_6.oncvpsp3.dojo.v0.4.1-std.upf",
"O.nc.z_6.oncvpsp3.dojo.v0.4.1-str.upf",
"O.nc.z_6.oncvpsp4.sg15.v0.upf",
"O.nc.z_6.oncvpsp4.spms.v1.upf",
"O.paw.z_6.atompaw.jth.v1.1-std.upf",
"O.paw.z_6.atompaw.jth.v1.1-str.upf",
"O.paw.z_6.ld1.psl.v0.1.upf",
"O.paw.z_6.ld1.psl.v1.0.0-high.upf",
"O.paw.z_6.ld1.psl.v1.0.0-low.upf",
"O.us.z_6.ld1.psl.v0.1.upf",
"O.us.z_6.ld1.psl.v1.0.0-high.upf",
"O.us.z_6.ld1.psl.v1.0.0-low.upf",
"O.us.z_6.uspp.gbrv.v1.2.upf",
]

#criteria = "efficiency"
criteria = "precision"
comment = "Nitrigen"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud-mem'
mpiprocs = 32 # 128 for eiger
npool = 4 # 8 for eiger
base_path = "./oxygen_psp/"

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment}  -- {pseudo_path}"
    os.system(command)
    # print(command)
    print(f"Launched {pseudo}")
