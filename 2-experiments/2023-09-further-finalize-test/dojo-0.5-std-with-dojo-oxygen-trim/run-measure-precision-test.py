import os

## Cutoff from recommendation of convergence test
#pseudos = {
# 'Bi.nc.z_15.oncvpsp4.dojo.v0.5.0-std-trim': (65.0, 165.0),
# 'Cd.nc.z_20.oncvpsp3.dojo.v0.5.0-std-trim': (120.0, 280.0),
# 'Cu.nc.z_19.oncvpsp3.dojo.v0.5.0-std-trim': (90.0, 280.0),
# 'Pb.nc.z_14.oncvpsp4.dojo.v0.5.0-std-trim': (120.0, 240.0),
# 'Pd.nc.z_18.oncvpsp3.dojo.v0.5.0-std-trim': (80.0, 200.0),
# 'Po.nc.z_16.oncvpsp4.dojo.v0.5.0-std-trim': (65.0, 175.0),
# 'Rn.nc.z_18.oncvpsp4.dojo.v0.5.0-std-trim': (120.0, 240.0),
# 'Te.nc.z_16.oncvpsp4.dojo.v0.5.0-std-trim': (75.0, 195.0),
# 'Tl.nc.z_13.oncvpsp4.dojo.v0.5.0-std-trim': (120.0, 240.0),
# 'Xe.nc.z_8.oncvpsp4.dojo.v0.5.0-std-trim': (150.0, 300.0),
# 'Zn.nc.z_20.oncvpsp3.dojo.v0.5.0-std-trim': (80.0, 280.0),
#}

# Trimed with exactly the same cutoff as dojo suggestion
pseudos = {
 'Na.nc.z_9.oncvpsp3.dojo.v0.5.0-std': (96.0, 96.0*4),
 #'Bi.nc.z_15.oncvpsp4.dojo.v0.5.0-std-trim': (74.0, 74.0*4),
 #'Cd.nc.z_20.oncvpsp3.dojo.v0.5.0-std-trim': (120.0, 280.0),
 #'Cu.nc.z_19.oncvpsp3.dojo.v0.5.0-std-trim': (90.0, 280.0),
 #'Pb.nc.z_14.oncvpsp4.dojo.v0.5.0-std-trim': (120.0, 240.0),
 #'Pd.nc.z_18.oncvpsp3.dojo.v0.5.0-std-trim': (80.0, 200.0),
 #'Po.nc.z_16.oncvpsp4.dojo.v0.5.0-std-trim': (76.0, 76.0*4),
 #'Rn.nc.z_18.oncvpsp4.dojo.v0.5.0-std-trim': (84.0, 84.0*4),
 #'Te.nc.z_16.oncvpsp4.dojo.v0.5.0-std-trim': (92.0, 92.0*4),
 #'Tl.nc.z_13.oncvpsp4.dojo.v0.5.0-std-trim': (74.0, 74.0*4),
 #'Xe.nc.z_8.oncvpsp4.dojo.v0.5.0-std-trim': (150.0, 300.0),
 #'Zn.nc.z_20.oncvpsp3.dojo.v0.5.0-std-trim': (80.0, 280.0),
}

comment = "dojo-0.5-std-trim"
computer = 'eiger-mc-mr33-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 8 for eiger
dojo_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.5-standard"
dojo_trim_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.5-standard-trim"
oxygen_pseudo_path = os.path.join(dojo_base_path, "O.nc.z_6.oncvpsp3.dojo.v0.5.0-std.upf")

#dojo_oxygen_ecutwfc = 80
#dojo_oxygen_ecutrho = 280

dojo_oxygen_ecutwfc = 96
dojo_oxygen_ecutrho = 96*4

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(dojo_trim_base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --oxygen-pseudo {oxygen_pseudo_path} --oxygen-ecutwfc {dojo_oxygen_ecutwfc} --oxygen-ecutrho {dojo_oxygen_ecutrho} --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    print(f"Launched {pseudo}")
    # print(command)