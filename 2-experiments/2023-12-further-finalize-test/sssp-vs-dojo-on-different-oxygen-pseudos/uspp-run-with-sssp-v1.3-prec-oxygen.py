import os

pseudos = {
 # GBRV pseudos
 'As.us.z_5.uspp.gbrv.v1': (60.0, 360.0),
 'Au.us.z_11.uspp.gbrv.v1': (65.0, 488.0),
 'Ca.us.z_10.uspp.gbrv.v1': (30.0, 180.0),
 'Ge.us.z_14.uspp.gbrv.v1.4': (40.0, 240.0),
 'Li.us.z_3.uspp.gbrv.v1.4': (40.0, 240.0),
 'Mg.us.z_10.uspp.gbrv.v1.4': (40.0, 240.0),
 'Mn.us.z_15.uspp.gbrv.v1.5': (150.0, 2400.0),
}

comment = "gbrv-with-sssp-v1.3-oxygen"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud-mem'
mpiprocs = 36 # 128 for eiger
npool = 4 # 8 for eiger
dojo_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.5-standard"
sssp_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/MIX-SSSP-precision-1.3.0-recollected"
uspp_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-GBRV-1.x"

oxygen_pseudo_path = os.path.join(sssp_base_path, 'O.paw.z_6.ld1.psl.v0.1.upf')

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    base_path = uspp_base_path
        
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)

    command = f"aiida-sssp-workflow launch --property measure.precision --oxygen-pseudo {oxygen_pseudo_path} --oxygen-ecutwfc 70 --oxygen-ecutrho 560 --configuration XO --configuration XO2 --configuration XO3 --configuration X2O --configuration X2O3 --configuration X2O5 --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    print(f"Launched {pseudo}")
    # print(command)
