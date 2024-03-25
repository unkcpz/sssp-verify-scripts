import os

pseudos = {
 #'Ag.paw.z_11.atompaw.jth.v1.1-std': (50.0, 400.0),
 #'Al.paw.z_3.atompaw.jth.v1.1-std': (35.0, 210.0),
 #'Ar.paw.z_8.atompaw.jth.v1.1-std': (75.0, 450.0),
 #'As.paw.z_15.atompaw.jth.v1.1-std': (35.0, 262.0),
 #'At.paw.z_7.atompaw.jth.v1.1-std': (30.0, 180.0),
 #'Au.paw.z_11.atompaw.jth.v1.1-std': (45.0, 270.0),
 #'B.paw.z_3.atompaw.jth.v1.1-std': (40.0, 240.0),
 #'Ba.paw.z_10.atompaw.jth.v1.1-std': (200.0, 1600.0),
 'Be.paw.z_4.atompaw.jth.v1.1-std': (65.0, 390.0),
 #'Bi.paw.z_15.atompaw.jth.v1.1-std': (60.0, 360.0),
 #'Br.paw.z_7.atompaw.jth.v1.1-std': (50.0, 375.0),
 #'C.paw.z_4.atompaw.jth.v1.1-std': (120.0, 900.0),
 #'Ca.paw.z_10.atompaw.jth.v1.1-std': (45.0, 270.0),
 #'Cd.paw.z_12.atompaw.jth.v1.1-std': (200.0, 1200.0),
 #'Cl.paw.z_7.atompaw.jth.v1.1-std': (40.0, 240.0),
 #'Co.paw.z_17.atompaw.jth.v1.1-std': (55.0, 990.0),
 #'Cr.paw.z_14.atompaw.jth.v1.1-std': (200.0, 3600.0),
 #'Cs.paw.z_9.atompaw.jth.v1.1-std': (200.0, 1600.0),
 #'Cu.paw.z_19.atompaw.jth.v1.1-std': (90.0, 720.0),
 #'F.paw.z_7.atompaw.jth.v1.1-std': (75.0, 450.0),
 #'Fe.paw.z_16.atompaw.jth.v1.1-std': (70.0, 1120.0),
 #'Ga.paw.z_13.atompaw.jth.v1.1-std': (30.0, 225.0),
 #'Ge.paw.z_14.atompaw.jth.v1.1-std': (50.0, 315.0),
 #'H.paw.z_1.atompaw.jth.v1.1-std': (45.0, 270.0),
 #'He.paw.z_2.atompaw.jth.v1.1-std': (120.0, 720.0),
 #'Hf.paw.z_12.atompaw.jth.v1.1-std': (50.0, 405.0),
 'Hg.paw.z_12.atompaw.jth.v1.1-std': (150.0, 1125.0),
 'I.paw.z_7.atompaw.jth.v1.1-std': (30.0, 180.0),
 #'In.paw.z_13.atompaw.jth.v1.1-std': (150.0, 900.0),
 'Ir.paw.z_15.atompaw.jth.v1.1-std': (35.0, 210.0),
 #'K.paw.z_9.atompaw.jth.v1.1-std': (50.0, 300.0),
 'Kr.paw.z_8.atompaw.jth.v1.1-std': (90.0, 675.0),
 #'Li.paw.z_3.atompaw.jth.v1.1-std': (90.0, 540.0),
 #'Mg.paw.z_10.atompaw.jth.v1.1-std': (120.0, 780.0),
 'Mn.paw.z_15.atompaw.jth.v1.1-std': (80.0, 640.0),
 'Mo.paw.z_14.atompaw.jth.v1.1-std': (55.0, 358.0),
 'N.paw.z_5.atompaw.jth.v1.1-std': (80.0, 640.0),
 #'Na.paw.z_9.atompaw.jth.v1.1-std': (120.0, 720.0),
 #'Nb.paw.z_13.atompaw.jth.v1.1-std': (45.0, 270.0),
 #'Ne.paw.z_8.atompaw.jth.v1.1-std': (200.0, 1600.0),
 #'Ni.paw.z_18.atompaw.jth.v1.1-std': (80.0, 640.0),
 #'O.paw.z_6.atompaw.jth.v1.1-std': (70.0, 560.0),
 #'Os.paw.z_8.atompaw.jth.v1.1-std': (30.0, 180.0),
 'P.paw.z_5.atompaw.jth.v1.1-std': (200.0, 1200.0),
 'Pb.paw.z_14.atompaw.jth.v1.1-std': (60.0, 360.0),
 #'Pd.paw.z_18.atompaw.jth.v1.1-std': (60.0, 480.0),
 #'Po.paw.z_6.atompaw.jth.v1.1-std': (30.0, 180.0),
 'Pt.paw.z_10.atompaw.jth.v1.1-std': (40.0, 240.0),
 #'Rb.paw.z_9.atompaw.jth.v1.1-std': (35.0, 210.0),
 #'Re.paw.z_15.atompaw.jth.v1.1-std': (65.0, 520.0),
 #'Rh.paw.z_17.atompaw.jth.v1.1-std': (120.0, 960.0),
 'Rn.paw.z_8.atompaw.jth.v1.1-std': (55.0, 412.0),
 #'Ru.paw.z_16.atompaw.jth.v1.1-std': (30.0, 210.0),
 #'S.paw.z_6.atompaw.jth.v1.1-std': (45.0, 270.0),
 'Sb.paw.z_15.atompaw.jth.v1.1-std': (40.0, 240.0),
 'Sc.paw.z_11.atompaw.jth.v1.1-std': (55.0, 440.0),
 #'Se.paw.z_6.atompaw.jth.v1.1-std': (30.0, 180.0),
 #'Si.paw.z_4.atompaw.jth.v1.1-std': (30.0, 240.0),
 #'Sn.paw.z_14.atompaw.jth.v1.1-std': (150.0, 900.0),
 #'Sr.paw.z_10.atompaw.jth.v1.1-std': (55.0, 440.0),
 'Ta.paw.z_13.atompaw.jth.v1.1-std': (65.0, 520.0),
 'Tc.paw.z_15.atompaw.jth.v1.1-std': (30.0, 240.0),
 'Te.paw.z_6.atompaw.jth.v1.1-std': (35.0, 210.0),
 #'Ti.paw.z_12.atompaw.jth.v1.1-std': (45.0, 338.0),
 #'Tl.paw.z_13.atompaw.jth.v1.1-std': (80.0, 480.0),
 'V.paw.z_13.atompaw.jth.v1.1-std': (50.0, 300.0),
 'W.paw.z_14.atompaw.jth.v1.1-std': (65.0, 520.0),
 #'Xe.paw.z_8.atompaw.jth.v1.1-std': (30.0, 240.0),
 'Y.paw.z_11.atompaw.jth.v1.1-std': (45.0, 360.0),
 'Zn.paw.z_12.atompaw.jth.v1.1-std': (60.0, 480.0),
 'Zr.paw.z_12.atompaw.jth.v1.1-std': (35.0, 210.0),
}

comment = "PAW-JTH-std-1.1"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-psi15-mem'
mpiprocs = 36 # 128 for eiger
npool = 4 # 8 for eiger
base_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-JTH1.1-standard'
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
