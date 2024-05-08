import os

pseudos = {
    #'N.paw.z_5.atompaw.jth.v1.1-std': (80.0, 600.0),
    #'N.nc.z_5.oncvpsp3.dojo.v0.4.1-std': (75.0, 188.0),
    'N.nc.z_5.oncvpsp4.sg15.v0': (75.0, 150.0),
    'N.nc.z_5.oncvpsp4.spms.v1': (75.0, 188.0),
    'N.paw.z_5.atompaw.jth.v1.1-str': (80.0, 600.0),
    'N.paw.z_5.ld1.psl.v0.1': (55.0, 330.0),
    'N.paw.z_5.ld1.psl.v1.0.0-high': (45.0, 270.0),
    'N.us.z_5.ld1.psl.v0.1': (55.0, 330.0),
    'N.us.z_5.ld1.psl.v1.0.0-high': (45.0, 270.0),
    'N.us.z_5.ld1.theose.v0': (50.0, 300.0),
    'N.us.z_5.uspp.gbrv.v1.2': (55.0, 330.0),
}

#criteria = "efficiency"
criteria = "precision"
comment = "Nitrigen"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud-mem'
mpiprocs = 32 # 128 for eiger
npool = 4 # 8 for eiger
base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/_sssp_pbe/N/"

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment}  -- {pseudo_path}"
    os.system(command)
    # print(command)
    print(f"Launched {pseudo}")
