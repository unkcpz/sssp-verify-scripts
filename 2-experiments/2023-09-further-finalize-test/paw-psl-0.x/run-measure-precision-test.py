import os

pseudos = {
 'Ag.paw.z_11.ld1.psl.v0.1': (200.0, 1600.0),
 'Al.paw.z_3.ld1.psl.v0.1': (40.0, 300.0),
 'Ar.paw.z_8.ld1.psl.v0.3.0': (90.0, 720.0),
 'As.paw.z_5.ld1.psl.v0.2': (35.0, 210.0),
 'Au.paw.z_11.ld1.psl.v0.3.0': (45.0, 360.0),
 'B.paw.z_3.ld1.psl.v0.1': (60.0, 360.0),
 'Bi.paw.z_15.ld1.psl.v0.2.2': (50.0, 300.0),
 'Br.paw.z_7.ld1.psl.v0.2': (35.0, 210.0),
 'C.paw.z_4.ld1.psl.v0.1': (60.0, 360.0),
 'Cd.paw.z_12.ld1.psl.v0.3.1': (90.0, 540.0),
 'Cl.paw.z_7.ld1.psl.v0.3.0': (40.0, 240.0),
 'Cu.paw.z_11.ld1.psl.v0.2': (150.0, 1050.0),
 'F.paw.z_7.ld1.psl.v0.1': (100.0, 750.0),
 'Fe.paw.z_16.ld1.psl.v0.2.1': (120.0, 960.0),
 'Ga.paw.z_13.ld1.psl.v0.2': (80.0, 640.0),
 'Ge.paw.z_14.ld1.psl.v0.3.1': (55.0, 400.0),
 'H.paw.z_1.ld1.psl.v0.1': (45.0, 270.0),
 'Hg.paw.z_12.ld1.psl.v0.2.2': (200.0, 1400.0),
 'In.paw.z_13.ld1.psl.v0.2.2': (50.0, 300.0),
 'Ir.paw.z_9.ld1.psl.v0.2.3': (150.0, 1200.0),
 'Li.paw.z_3.ld1.psl.v0.2.1': (80.0, 480.0),
 'Mo.paw.z_14.ld1.psl.v0.3.0': (60.0, 480.0),
 'N.paw.z_5.ld1.psl.v0.1': (55.0, 330.0),
 'Na.paw.z_9.ld1.psl.v0.2': (60.0, 360.0),
 'Nb.paw.z_13.ld1.psl.v0.3.0': (50.0, 300.0),
 'Ni.paw.z_10.ld1.psl.v0.1': (55.0, 450.0),
 'O.paw.z_6.ld1.psl.v0.1': (70.0, 560.0),
 'P.paw.z_5.ld1.psl.v0.1': (35.0, 210.0),
 'Pb.paw.z_14.ld1.psl.v0.2.2': (50.0, 300.0),
 'Pd.paw.z_10.ld1.psl.v0.3.0': (90.0, 675.0),
 'Pt.paw.z_10.ld1.psl.v0.1': (55.0, 330.0),
 'Rh.paw.z_17.ld1.psl.v0.3.0': (75.0, 480.0),
 'S.paw.z_6.ld1.psl.v0.1': (40.0, 240.0),
 'Sc.paw.z_11.ld1.psl.v0.2.3': (75.0, 450.0),
 'Se.paw.z_6.ld1.psl.v0.2': (35.0, 210.0),
 'Si.paw.z_4.ld1.psl.v0.1': (30.0, 180.0),
 'Sr.paw.z_10.ld1.psl.v0.2.3': (50.0, 400.0),
 'Ta.paw.z_13.ld1.psl.v0.2': (70.0, 560.0),
 'Tc.paw.z_15.ld1.psl.v0.3.0': (75.0, 600.0),
 'Tl.paw.z_13.ld1.psl.v0.2.3': (65.0, 390.0),
 'Zr.paw.z_12.ld1.psl.v0.2.3': (50.0, 300.0),
}

comment = "PAW-PSL-0.x"
computer = 'eiger-mc-mr32-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 16 for eiger
base_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-PSL0.x'
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
