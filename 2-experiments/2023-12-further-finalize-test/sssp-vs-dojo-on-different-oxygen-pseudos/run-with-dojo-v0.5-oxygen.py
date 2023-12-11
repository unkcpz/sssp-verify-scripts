import os

pseudos = {
 # DOJO pseudos
 #'As.nc.z_15.oncvpsp3.dojo.v0.5.0-std': (80.0, 245.0),
 'Au.nc.z_19.oncvpsp3.dojo.v0.5.0-std': (70.0, 175.0),
 'Ca.nc.z_10.oncvpsp3.dojo.v0.5.0-std': (70.0, 175.0),
 'Li.nc.z_3.oncvpsp3.dojo.v0.5.0-std': (75.0, 188.0),
 'Mg.nc.z_10.oncvpsp3.dojo.v0.5.0-std': (120.0, 480.0),
 'Ge.nc.z_14.oncvpsp3.dojo.v0.5.0-std': (75.0, 228.0),
 'Mn.nc.z_15.oncvpsp3.dojo.v0.5.0-std': (150.0, 300.0),

 # SSSP pseudos
 'Br.us.z_7.uspp.gbrv.v1.4': (30.0, 180.0),
 'Te.us.z_6.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Fe.paw.z_16.ld1.psl.v0.2.1': (120.0, 960.0),
 'Hg.us.z_12.uspp.gbrv.v1': (150.0, 900.0),
}

comment = "with-dojo-v0.5-oxygen"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud-mem'
mpiprocs = 36 # 128 for eiger
npool = 4 # 8 for eiger
dojo_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.5-standard"
sssp_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/MIX-SSSP-precision-1.3.0-recollected"

oxygen_pseudo_path = os.path.join(dojo_base_path, 'O.nc.z_6.oncvpsp3.dojo.v0.5.0-std.upf')

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    if "dojo" in pseudo:
        base_path = dojo_base_path
    else:
        base_path = sssp_base_path
        
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)

    command = f"aiida-sssp-workflow launch --property measure.precision --oxygen-pseudo {oxygen_pseudo_path} --oxygen-ecutwfc 80 --oxygen-ecutrho 280 --configuration XO --configuration XO2 --configuration XO3 --configuration X2O --configuration X2O3 --configuration X2O5 --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    print(f"Launched {pseudo}")
    # print(command)
