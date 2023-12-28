import os

pseudos = {
 'Ag.nc.z_19.oncvpsp4.sg15.v0': (50.0, 140.0),
 'Al.nc.z_11.oncvpsp4.sg15.v0': (70.0, 140.0),
 'As.nc.z_5.oncvpsp4.sg15.v0': (30.0, 60.0),
 'Au.nc.z_19.oncvpsp4.sg15.v0': (65.0, 130.0),
 'B.nc.z_3.oncvpsp4.sg15.v0': (40.0, 88.0),
 'Ba.nc.z_10.oncvpsp4.sg15.v0': (45.0, 90.0),
 'Be.nc.z_4.oncvpsp4.sg15.v0': (50.0, 135.0),
 'Bi.nc.z_15.oncvpsp4.sg15.v0': (40.0, 120.0),
 'Br.nc.z_7.oncvpsp4.sg15.v0': (200.0, 800.0),
 'C.nc.z_4.oncvpsp4.sg15.v0': (75.0, 150.0),
 'Ca.nc.z_10.oncvpsp4.sg15.v0': (40.0, 120.0),
 'Cd.nc.z_20.oncvpsp4.sg15.v0': (40.0, 88.0),
 'Cl.nc.z_7.oncvpsp4.sg15.v0': (50.0, 100.0),
 'Co.nc.z_17.oncvpsp4.sg15.v0': (60.0, 220.0),
 'Cr.nc.z_14.oncvpsp4.sg15.v0': (50.0, 125.0),
 'Cs.nc.z_9.oncvpsp4.sg15.v0': (45.0, 90.0),
 'Cu.nc.z_19.oncvpsp4.sg15.v0': (80.0, 245.0),
 'F.nc.z_7.oncvpsp4.sg15.v0': (150.0, 600.0),
 'Fe.nc.z_16.oncvpsp4.sg15.v0': (70.0, 228.0),
 'Ga.nc.z_13.oncvpsp4.sg15.v0': (150.0, 300.0),
 'Ge.nc.z_14.oncvpsp4.sg15.v0': (45.0, 180.0),
 'H.nc.z_1.oncvpsp4.sg15.v0': (70.0, 140.0),
 'He.nc.z_2.oncvpsp4.sg15.v0': (80.0, 320.0),
 'Hf.nc.z_26.oncvpsp4.sg15.v0': (65.0, 158.0),
 'Hg.nc.z_20.oncvpsp4.sg15.v0': (200.0, 800.0),
 'I.nc.z_17.oncvpsp4.sg15.v0': (65.0, 160.0),
 'In.nc.z_13.oncvpsp4.sg15.v0': (70.0, 140.0),
 'Ir.nc.z_17.oncvpsp4.sg15.v0': (200.0, 800.0),
 'K.nc.z_9.oncvpsp4.sg15.v0': (65.0, 260.0),
 'Kr.nc.z_8.oncvpsp4.sg15.v0': (60.0, 240.0),
 'Li.nc.z_3.oncvpsp4.sg15.v0': (65.0, 150.0),
 'Mg.nc.z_10.oncvpsp4.sg15.v0': (100.0, 200.0),
 'Mn.nc.z_15.oncvpsp4.sg15.v0': (65.0, 150.0),
 'Mo.nc.z_14.oncvpsp4.sg15.v0': (45.0, 90.0),
 'N.nc.z_5.oncvpsp4.sg15.v0': (75.0, 210.0),
 'Na.nc.z_9.oncvpsp4.sg15.v0': (100.0, 200.0),
 'Nb.nc.z_13.oncvpsp4.sg15.v0': (60.0, 135.0),
 'Ne.nc.z_8.oncvpsp4.sg15.v0': (65.0, 260.0),
 'Ni.nc.z_18.oncvpsp4.sg15.v0': (200.0, 700.0),
 'O.nc.z_6.oncvpsp4.sg15.v0': (100.0, 300.0),
 'Os.nc.z_16.oncvpsp4.sg15.v0': (60.0, 200.0),
 'P.nc.z_5.oncvpsp4.sg15.v0': (65.0, 130.0),
 'Pb.nc.z_14.oncvpsp4.sg15.v0': (35.0, 120.0),
 'Pd.nc.z_18.oncvpsp4.sg15.v0': (50.0, 140.0),
 'Pt.nc.z_18.oncvpsp4.sg15.v0': (45.0, 120.0),
 'Rb.nc.z_9.oncvpsp4.sg15.v0': (30.0, 60.0),
 'Re.nc.z_15.oncvpsp4.sg15.v0': (45.0, 100.0),
 'Rh.nc.z_17.oncvpsp4.sg15.v0': (60.0, 120.0),
 'Ru.nc.z_16.oncvpsp4.sg15.v0': (45.0, 112.0),
 'S.nc.z_6.oncvpsp4.sg15.v0': (75.0, 150.0),
 'Sb.nc.z_15.oncvpsp4.sg15.v0': (55.0, 160.0),
 'Sc.nc.z_11.oncvpsp4.sg15.v0': (45.0, 90.0),
 'Se.nc.z_6.oncvpsp4.sg15.v0': (30.0, 60.0),
 'Si.nc.z_4.oncvpsp4.sg15.v0': (45.0, 90.0),
 'Sn.nc.z_14.oncvpsp4.sg15.v0': (45.0, 140.0),
 'Sr.nc.z_10.oncvpsp4.sg15.v0': (30.0, 60.0),
 'Ta.nc.z_27.oncvpsp4.sg15.v0': (50.0, 158.0),
 'Tc.nc.z_15.oncvpsp4.sg15.v0': (45.0, 112.0),
 'Te.nc.z_16.oncvpsp4.sg15.v0': (55.0, 180.0),
 'Ti.nc.z_12.oncvpsp4.sg15.v0': (50.0, 135.0),
 'Tl.nc.z_13.oncvpsp4.sg15.v0': (35.0, 120.0),
 'V.nc.z_13.oncvpsp4.sg15.v0': (200.0, 400.0),
 'W.nc.z_28.oncvpsp4.sg15.v0': (100.0, 400.0),
 'Xe.nc.z_18.oncvpsp4.sg15.v0': (75.0, 300.0),
 'Y.nc.z_11.oncvpsp4.sg15.v0': (45.0, 90.0),
 'Zn.nc.z_20.oncvpsp4.sg15.v0': (80.0, 210.0),
 'Zr.nc.z_12.oncvpsp4.sg15.v0': (50.0, 100.0),
}


comment = "sg15-with-dojo-oxygen"
computer = 'eiger-mc-mr33-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 8 for eiger
base_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-SG15-ONCVPSP4'
dojo_base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-DOJOv0.5-standard"
oxygen_pseudo_path = os.path.join(dojo_base_path, "O.nc.z_6.oncvpsp3.dojo.v0.5.0-std.upf")

dojo_oxygen_ecutwfc = 80
dojo_oxygen_ecutrho = 280

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --oxygen-pseudo {oxygen_pseudo_path} --oxygen-ecutwfc {dojo_oxygen_ecutwfc} --oxygen-ecutrho {dojo_oxygen_ecutrho} --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    print(f"Launched {pseudo}")
    # print(command)
