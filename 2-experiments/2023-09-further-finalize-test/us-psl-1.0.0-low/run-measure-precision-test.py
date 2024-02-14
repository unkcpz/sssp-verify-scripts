import os

pseudos = {'Ag.us.z_11.ld1.psl.v1.0.0-low': (70.0, 420.0),
 'Al.us.z_3.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Ar.us.z_8.ld1.psl.v1.0.0-low': (120.0, 720.0),
 'As.us.z_5.ld1.psl.v1.0.0-low': (35.0, 210.0),
 'Au.us.z_11.ld1.psl.v1.0.0-low.n': (45.0, 338.0),
 'Au.us.z_19.ld1.psl.v1.0.0-low.spn': (90.0, 540.0),
 'Be.us.z_2.ld1.psl.v1.0.0-low.n': (30.0, 180.0),
 'Be.us.z_4.ld1.psl.v1.0.0-low.sl': (150.0, 1050.0),
 'Br.us.z_7.ld1.psl.v1.0.0-low': (40.0, 240.0),
 'Cd.us.z_12.ld1.psl.v1.0.0-low': (90.0, 640.0),
 'Cl.us.z_7.ld1.psl.v1.0.0-low': (65.0, 390.0),
 'Co.us.z_9.ld1.psl.v1.0.0-low': (75.0, 600.0),
 'Cs.us.z_9.ld1.psl.v1.0.0-low': (50.0, 300.0),
 'Cu.us.z_11.ld1.psl.v1.0.0-low': (70.0, 455.0),
 'Fe.us.z_8.ld1.psl.v1.0.0-low': (75.0, 600.0),
 'Ga.us.z_13.ld1.psl.v1.0.0-low': (65.0, 390.0),
 'Ge.us.z_4.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Hf.us.z_12.ld1.psl.v1.0.0-low': (60.0, 480.0),
 'Hg.us.z_12.ld1.psl.v1.0.0-low': (120.0, 720.0),
 'I.us.z_7.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Ir.us.z_17.ld1.psl.v1.0.0-low.spn': (70.0, 520.0),
 'Ir.us.z_9.ld1.psl.v1.0.0-low.n': (45.0, 270.0),
 'Li.us.z_3.ld1.psl.v1.0.0-low': (75.0, 450.0),
 'Mg.us.z_10.ld1.psl.v1.0.0-low': (150.0, 900.0),
 'Na.us.z_9.ld1.psl.v1.0.0-low': (150.0, 900.0),
 'Ni.us.z_10.ld1.psl.v1.0.0-low': (55.0, 450.0),
 'O.us.z_6.ld1.psl.v1.0.0-low': (65.0, 520.0),
 'Os.us.z_16.ld1.psl.v1.0.0-low': (60.0, 360.0),
 'P.us.z_5.ld1.psl.v1.0.0-low': (45.0, 270.0),
 'Pd.us.z_10.ld1.psl.v1.0.0-low': (55.0, 440.0),
 'Pt.us.z_10.ld1.psl.v1.0.0-low.n': (55.0, 330.0),
 'Pt.us.z_18.ld1.psl.v1.0.0-low.spn': (80.0, 560.0),
 'Re.us.z_15.ld1.psl.v1.0.0-low': (60.0, 360.0),
 'S.us.z_6.ld1.psl.v1.0.0-low': (50.0, 300.0),
 'Sb.us.z_5.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Se.us.z_6.ld1.psl.v1.0.0-low': (35.0, 210.0),
 'Si.us.z_4.ld1.psl.v1.0.0-low': (40.0, 240.0),
 'Ta.us.z_13.ld1.psl.v1.0.0-low': (50.0, 375.0),
 'Te.us.z_6.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'V.us.z_13.ld1.psl.v1.0.0-low': (150.0, 900.0),
 'W.us.z_14.ld1.psl.v1.0.0-low': (120.0, 900.0),
 'Zn.us.z_12.ld1.psl.v1.0.0-low': (90.0, 720.0)}

comment = "US-PSL-1.0.0-low"
computer = 'eiger-mc-mr32-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 16 for eiger
base_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-PSL1.0.0-low'
dojo_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.5-standard"
sssp_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/MIX-SSSP-precision-1.3.0-recollected"
#oxygen_pseudo_path = os.path.join(dojo_base_path, "O.nc.z_6.oncvpsp3.dojo.v0.5.0-std.upf")
oxygen_pseudo_path = os.path.join(sssp_base_path, 'O.paw.z_6.ld1.psl.v0.1.upf')

oxygen_ecutwfc = 70
oxygen_ecutrho = 560

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --oxygen-pseudo {oxygen_pseudo_path} --oxygen-ecutwfc {oxygen_ecutwfc} --oxygen-ecutrho {oxygen_ecutrho} --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    print(f"Launched {pseudo}")
    # print(command)
