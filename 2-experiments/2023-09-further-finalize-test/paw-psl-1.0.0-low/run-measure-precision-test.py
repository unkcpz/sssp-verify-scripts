import os

pseudos = {
 'Ag.paw.z_11.ld1.psl.v1.0.0-low': (70.0, 420.0),
 'Al.paw.z_3.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Ar.paw.z_8.ld1.psl.v1.0.0-low': (120.0, 720.0),
 'As.paw.z_5.ld1.psl.v1.0.0-low': (35.0, 210.0),
 'Au.paw.z_11.ld1.psl.v1.0.0-low.n': (45.0, 338.0),
 'Au.paw.z_19.ld1.psl.v1.0.0-low.spn': (90.0, 540.0),
 'Be.paw.z_2.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Be.paw.z_4.ld1.psl.v1.0.0-low.sl': (150.0, 1050.0),
 'Br.paw.z_7.ld1.psl.v1.0.0-low': (40.0, 240.0),
 'Cd.paw.z_12.ld1.psl.v1.0.0-low': (90.0, 640.0),
 'Cl.paw.z_7.ld1.psl.v1.0.0-low': (65.0, 390.0),
 'Co.paw.z_9.ld1.psl.v1.0.0-low': (75.0, 720.0),
 'Cs.paw.z_9.ld1.psl.v1.0.0-low': (50.0, 300.0),
 'Cu.paw.z_11.ld1.psl.v1.0.0-low': (70.0, 455.0),
 'Fe.paw.z_8.ld1.psl.v1.0.0-low': (75.0, 600.0),
 'Ga.paw.z_13.ld1.psl.v1.0.0-low': (65.0, 390.0),
 'Ge.paw.z_4.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Hf.paw.z_12.ld1.psl.v1.0.0-low': (60.0, 480.0),
 'Hg.paw.z_12.ld1.psl.v1.0.0-low': (120.0, 720.0),
 'I.paw.z_7.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Ir.paw.z_17.ld1.psl.v1.0.0-low.spn': (70.0, 520.0),
 'Ir.paw.z_9.ld1.psl.v1.0.0-low.n': (45.0, 270.0),
 'Li.paw.z_3.ld1.psl.v1.0.0-low': (75.0, 450.0),
 'Mg.paw.z_10.ld1.psl.v1.0.0-low': (150.0, 900.0),
 'Na.paw.z_9.ld1.psl.v1.0.0-low': (150.0, 900.0),
 'Ni.paw.z_10.ld1.psl.v1.0.0-low': (55.0, 450.0),
 'O.paw.z_6.ld1.psl.v1.0.0-low': (65.0, 520.0),
 'Os.paw.z_16.ld1.psl.v1.0.0-low': (60.0, 360.0),
 'P.paw.z_5.ld1.psl.v1.0.0-low': (45.0, 270.0),
 'Pd.paw.z_10.ld1.psl.v1.0.0-low': (55.0, 440.0),
 'Pt.paw.z_10.ld1.psl.v1.0.0-low.n': (55.0, 330.0),
 'Pt.paw.z_18.ld1.psl.v1.0.0-low.spn': (80.0, 560.0),
 'Re.paw.z_15.ld1.psl.v1.0.0-low': (60.0, 360.0),
 'S.paw.z_6.ld1.psl.v1.0.0-low': (50.0, 300.0),
 'Sb.paw.z_5.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'Se.paw.z_6.ld1.psl.v1.0.0-low': (35.0, 210.0),
 'Si.paw.z_4.ld1.psl.v1.0.0-low': (40.0, 240.0),
 'Ta.paw.z_13.ld1.psl.v1.0.0-low': (50.0, 375.0),
 'Te.paw.z_6.ld1.psl.v1.0.0-low': (30.0, 180.0),
 'V.paw.z_13.ld1.psl.v1.0.0-low': (80.0, 490.0),
 'W.paw.z_14.ld1.psl.v1.0.0-low': (55.0, 440.0),
 'Zn.paw.z_12.ld1.psl.v1.0.0-low': (65.0, 390.0),
}

comment = "PAW-PSL-1.0.0-low"
computer = 'eiger-mc-mr33-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 16 for eiger
base_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-PSL1.0.0-low'
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
