import os

pseudos = {
 #'Ag.nc.z_19.oncvpsp3.dojo.v0.5.0-std': (90.0, 225.0),
 #'Al.nc.z_3.oncvpsp3.dojo.v0.5.0-std': (30.0, 60.0),
 #'Ar.nc.z_8.oncvpsp3.dojo.v0.5.0-std': (70.0, 165.0),
 #'As.nc.z_15.oncvpsp3.dojo.v0.5.0-std': (80.0, 245.0),
 #'Au.nc.z_19.oncvpsp3.dojo.v0.5.0-std': (70.0, 175.0),
 #'B.nc.z_3.oncvpsp3.dojo.v0.5.0-std': (65.0, 130.0),
 #'Ba.nc.z_10.oncvpsp4.dojo.v0.5.0-std': (55.0, 110.0),
 #'Be.nc.z_4.oncvpsp3.dojo.v0.5.0-std': (90.0, 240.0),
 #'Bi.nc.z_15.oncvpsp4.dojo.v0.5.0-std': (65.0, 165.0),
 #'Br.nc.z_7.oncvpsp3.dojo.v0.5.0-std': (45.0, 90.0),
 #'C.nc.z_4.oncvpsp3.dojo.v0.5.0-std': (75.0, 150.0),
 #'Ca.nc.z_10.oncvpsp3.dojo.v0.5.0-std': (70.0, 175.0),
 #'Cd.nc.z_20.oncvpsp3.dojo.v0.5.0-std': (120.0, 280.0),
 #'Cl.nc.z_7.oncvpsp3.dojo.v0.5.0-std': (55.0, 110.0),
 #'Co.nc.z_17.oncvpsp3.dojo.v0.5.0-std': (120.0, 240.0),
 #'Cr.nc.z_14.oncvpsp3.dojo.v0.5.0-std': (90.0, 270.0),
 #'Cs.nc.z_9.oncvpsp3.dojo.v0.5.0-std': (60.0, 120.0),
 #'Cu.nc.z_19.oncvpsp3.dojo.v0.5.0-std': (90.0, 280.0),
 #'F.nc.z_7.oncvpsp3.dojo.v0.5.0-std': (80.0, 260.0),
 #'Fe.nc.z_16.oncvpsp3.dojo.v0.5.0-std': (90.0, 270.0),
 #'Ga.nc.z_13.oncvpsp3.dojo.v0.5.0-std': (80.0, 245.0),
 #'Ge.nc.z_14.oncvpsp3.dojo.v0.5.0-std': (75.0, 228.0),
 #'H.nc.z_1.oncvpsp3.dojo.v0.5.0-std': (65.0, 130.0),
 #'He.nc.z_2.oncvpsp3.dojo.v0.5.0-std': (90.0, 320.0),
 #'Hf.nc.z_12.oncvpsp3.dojo.v0.5.0-std': (55.0, 110.0),
 #'Hg.nc.z_20.oncvpsp3.dojo.v0.5.0-std': (70.0, 175.0),
 #'I.nc.z_7.oncvpsp4.dojo.v0.5.0-std': (65.0, 130.0),
 #'In.nc.z_13.oncvpsp3.dojo.v0.5.0-std': (70.0, 180.0),
 #'Ir.nc.z_17.oncvpsp3.dojo.v0.5.0-std': (60.0, 150.0),
 #'K.nc.z_9.oncvpsp3.dojo.v0.5.0-std': (75.0, 158.0),
 #'Kr.nc.z_8.oncvpsp3.dojo.v0.5.0-std': (120.0, 240.0),
 #'Li.nc.z_3.oncvpsp3.dojo.v0.5.0-std': (75.0, 188.0),
 #'Mg.nc.z_10.oncvpsp3.dojo.v0.5.0-std': (120.0, 480.0),
 #'Mn.nc.z_15.oncvpsp3.dojo.v0.5.0-std': (150.0, 300.0),
 #'Mo.nc.z_14.oncvpsp3.dojo.v0.5.0-std': (75.0, 150.0),
 #'N.nc.z_5.oncvpsp3.dojo.v0.5.0-std': (75.0, 188.0),
 'Na.nc.z_9.oncvpsp3.dojo.v0.5.0-std': (100.0, 400.0),
 #'Nb.nc.z_13.oncvpsp3.dojo.v0.5.0-std': (75.0, 180.0),
 #'Ne.nc.z_8.oncvpsp3.dojo.v0.5.0-std': (200.0, 800.0),
 #'Ni.nc.z_18.oncvpsp3.dojo.v0.5.0-std': (200.0, 400.0),
 #'O.nc.z_6.oncvpsp3.dojo.v0.5.0-std': (80.0, 280.0),
 #'Os.nc.z_16.oncvpsp3.dojo.v0.5.0-std': (65.0, 150.0),
 #'P.nc.z_5.oncvpsp3.dojo.v0.5.0-std': (40.0, 80.0),
 #'Pb.nc.z_14.oncvpsp4.dojo.v0.5.0-std': (120.0, 240.0),
 #'Pd.nc.z_18.oncvpsp3.dojo.v0.5.0-std': (80.0, 200.0),
 #'Po.nc.z_16.oncvpsp4.dojo.v0.5.0-std': (65.0, 175.0),
 #'Pt.nc.z_18.oncvpsp3.dojo.v0.5.0-std': (80.0, 188.0),
 #'Rb.nc.z_9.oncvpsp4.dojo.v0.5.0-std': (50.0, 100.0),
 #'Re.nc.z_15.oncvpsp3.dojo.v0.5.0-std': (90.0, 180.0),
 #'Rh.nc.z_17.oncvpsp3.dojo.v0.5.0-std': (90.0, 192.0),
 #'Rn.nc.z_18.oncvpsp4.dojo.v0.5.0-std': (120.0, 240.0),
 #'Ru.nc.z_16.oncvpsp3.dojo.v0.5.0-std': (75.0, 175.0),
 #'S.nc.z_6.oncvpsp4.dojo.v0.5.0-std': (45.0, 90.0),
 #'Sb.nc.z_15.oncvpsp3.dojo.v0.5.0-std': (80.0, 228.0),
 #'Sc.nc.z_11.oncvpsp3.dojo.v0.5.0-std': (80.0, 200.0),
 #'Se.nc.z_16.oncvpsp3.dojo.v0.5.0-std': (80.0, 245.0),
 #'Si.nc.z_4.oncvpsp3.dojo.v0.5.0-std': (30.0, 60.0),
 #'Sn.nc.z_14.oncvpsp3.dojo.v0.5.0-std': (70.0, 210.0),
 #'Sr.nc.z_10.oncvpsp3.dojo.v0.5.0-std': (65.0, 130.0),
 #'Ta.nc.z_13.oncvpsp3.dojo.v0.5.0-std': (50.0, 100.0),
 #'Tc.nc.z_15.oncvpsp3.dojo.v0.5.0-std': (75.0, 175.0),
 #'Te.nc.z_16.oncvpsp4.dojo.v0.5.0-std': (75.0, 195.0),
 #'Ti.nc.z_12.oncvpsp3.dojo.v0.5.0-std': (80.0, 210.0),
 #'Tl.nc.z_13.oncvpsp4.dojo.v0.5.0-std': (120.0, 240.0),
 #'V.nc.z_13.oncvpsp3.dojo.v0.5.0-std': (80.0, 262.0),
 #'W.nc.z_14.oncvpsp3.dojo.v0.5.0-std': (100.0, 200.0),
 #'Xe.nc.z_8.oncvpsp4.dojo.v0.5.0-std': (150.0, 300.0),
 #'Y.nc.z_11.oncvpsp3.dojo.v0.5.0-std': (70.0, 140.0),
 #'Zn.nc.z_20.oncvpsp3.dojo.v0.5.0-std': (80.0, 280.0),
 #'Zr.nc.z_12.oncvpsp3.dojo.v0.5.0-std': (60.0, 135.0),
}

comment = "Na-cutoff-check"
computer = 'eiger-mc-mr33-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 16 for eiger
dojo_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.5-standard"
oxygen_pseudo_path = os.path.join(dojo_base_path, "O.nc.z_6.oncvpsp3.dojo.v0.5.0-std.upf")

dojo_oxygen_ecutwfc = 80
dojo_oxygen_ecutrho = 280

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(dojo_base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --configuration Diamond --oxygen-pseudo {oxygen_pseudo_path} --oxygen-ecutwfc {dojo_oxygen_ecutwfc} --oxygen-ecutrho {dojo_oxygen_ecutrho} --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    print(f"Launched {pseudo}")
    # print(command)