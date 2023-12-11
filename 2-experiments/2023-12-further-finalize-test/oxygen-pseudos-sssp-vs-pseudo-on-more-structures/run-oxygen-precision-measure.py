import os

# Get by runnning `verdi run ../../utils/get_recommened_cutoffs.py element-Oxygen/convergence/precision`
pseudos = {
 #'O.nc.z_6.oncvpsp3.dojo.v0.4.1-std': (80.0, 280.0),
 'O.nc.z_6.oncvpsp3.dojo.v0.4.1-str': (90.0, 300.0),
 'O.nc.z_6.oncvpsp4.sg15.v0': (100.0, 300.0),
 'O.nc.z_6.oncvpsp4.spms.v1': (80.0, 300.0),
 'O.paw.z_6.atompaw.jth.v1.1-std': (70.0, 560.0),
 'O.paw.z_6.atompaw.jth.v1.1-str': (70.0, 560.0),
 #'O.paw.z_6.ld1.psl.v0.1': (70.0, 560.0),
 'O.paw.z_6.ld1.psl.v1.0.0-high': (75.0, 600.0),
 'O.paw.z_6.ld1.psl.v1.0.0-low': (65.0, 520.0),
 #'O.us.z_7.ld1.psl.v0.1': (70.0, 990.0),
 'O.us.z_6.ld1.psl.v1.0.0-high': (75.0, 600.0),
 'O.us.z_6.ld1.psl.v1.0.0-low': (65.0, 520.0),
}


comment = "Oxygen"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud-mem'
mpiprocs = 36 # 128 for eiger
npool = 4 # 8 for eiger
base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/_sssp_pbe/O/"

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    # print(command)
