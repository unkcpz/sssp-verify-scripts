import os

pseudos = {
 'Ag.us.z_11.ld1.psl.v0.1': (200.0, 1600.0),
 'Al.us.z_3.ld1.psl.v0.1': (40.0, 300.0),
 'Ar.us.z_8.ld1.psl.v0.3.0': (120.0, 900.0),
 'As.us.z_5.ld1.psl.v0.2': (35.0, 210.0),
 'Au.us.z_11.ld1.psl.v0.3.0': (45.0, 360.0),
 'B.us.z_3.ld1.psl.v0.1': (70.0, 560.0),
 'Bi.us.z_15.ld1.psl.v0.2.2': (70.0, 560.0),
 'Br.us.z_7.ld1.psl.v0.2': (35.0, 210.0),
 'C.us.z_4.ld1.psl.v0.1': (60.0, 360.0),
 'Cd.us.z_12.ld1.psl.v0.3.1': (90.0, 540.0),
 'Cl.us.z_7.ld1.psl.v0.3.0': (40.0, 240.0),
 'Cu.us.z_11.ld1.psl.v0.2': (100.0, 800.0),
 'F.us.z_7.ld1.psl.v0.1': (100.0, 750.0),
 'Fe.us.z_16.ld1.psl.v0.2.1': (120.0, 960.0),
 'Ga.us.z_13.ld1.psl.v0.2': (50.0, 375.0),
 'Ge.us.z_14.ld1.psl.v0.3.1': (55.0, 400.0),
 'H.us.z_1.ld1.psl.v0.1': (45.0, 270.0),
 'Hg.us.z_12.ld1.psl.v0.2.2': (120.0, 720.0),
 'In.us.z_13.ld1.psl.v0.2.2': (50.0, 300.0),
 'Ir.us.z_9.ld1.psl.v0.2.3': (45.0, 270.0),
 'Li.us.z_3.ld1.psl.v0.2.1': (80.0, 480.0),
 'Mo.us.z_14.ld1.psl.v0.3.0': (60.0, 480.0),
 'N.us.z_5.ld1.psl.v0.1': (55.0, 330.0),
 'Na.us.z_9.ld1.psl.v0.2': (60.0, 360.0),
 'Nb.us.z_13.ld1.psl.v0.3.0': (50.0, 300.0),
 'Ni.us.z_10.ld1.psl.v0.1': (55.0, 450.0),
 'O.us.z_6.ld1.psl.v0.1': (70.0, 990.0),
 'P.us.z_5.ld1.psl.v0.1': (35.0, 210.0),
 'Pb.us.z_14.ld1.psl.v0.2.2': (45.0, 270.0),
 'Pd.us.z_10.ld1.psl.v0.3.0': (90.0, 675.0),
 'Pt.us.z_10.ld1.psl.v0.1': (55.0, 330.0),
 'Rh.us.z_17.ld1.psl.v0.3.0': (75.0, 480.0),
 'S.us.z_6.ld1.psl.v0.1': (40.0, 240.0),
 'Sc.us.z_11.ld1.psl.v0.2.3': (75.0, 450.0),
 'Se.us.z_6.ld1.psl.v0.2': (35.0, 210.0),
 'Si.us.z_4.ld1.psl.v0.1': (30.0, 180.0),
 'Sr.us.z_10.ld1.psl.v0.2.3': (50.0, 400.0),
 'Ta.us.z_13.ld1.psl.v0.2': (70.0, 560.0),
 'Tc.us.z_15.ld1.psl.v0.3.0': (75.0, 600.0),
 'Tl.us.z_13.ld1.psl.v0.2.3': (65.0, 390.0),
 'Zr.us.z_12.ld1.psl.v0.2.3': (50.0, 300.0),
}

comment = "us-psl-0.x"
computer = 'eiger-mc-mr32-mem'
#computer = 'daint-mc-mrcloud-mem'
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 8 # 8 for eiger
base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-PSL0.x"

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    print(f"Launched {pseudo}")
    # print(command)
