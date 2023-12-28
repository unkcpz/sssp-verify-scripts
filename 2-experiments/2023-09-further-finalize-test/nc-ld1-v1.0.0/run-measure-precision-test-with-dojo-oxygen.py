import os

pseudos = {
 'Ag.nc.z_11.ld1.psl.v1.0.0': (200.0, 800.0),
 'Al.nc.z_3.ld1.psl.v1.0.0': (60.0, 240.0),
 'Ar.nc.z_8.ld1.psl.v1.0.0': (200.0, 800.0),
 'As.nc.z_5.ld1.psl.v1.0.0': (200.0, 700.0),
 'At.nc.z_7.ld1.psl.v1.0.0': (150.0, 600.0),
 'Au.nc.z_11.ld1.psl.v1.0.0': (200.0, 800.0),
 'B.nc.z_3.ld1.psl.v1.0.0': (55.0, 220.0),
 'Ba.nc.z_2.ld1.psl.v1.0.0': (200.0, 800.0),
 'Be.nc.z_2.ld1.psl.v1.0.0': (60.0, 240.0),
 'Bi.nc.z_5.ld1.psl.v1.0.0': (90.0, 360.0),
 'Br.nc.z_7.ld1.psl.v1.0.0': (200.0, 800.0),
 'C.nc.z_4.ld1.psl.v1.0.0': (120.0, 240.0),
 'Ca.nc.z_2.ld1.psl.v1.0.0': (30.0, 75.0),
 'Cd.nc.z_12.ld1.psl.v1.0.0': (200.0, 800.0),
 'Cl.nc.z_7.ld1.psl.v1.0.0': (200.0, 800.0),
 'Co.nc.z_9.ld1.psl.v1.0.0': (200.0, 800.0),
 'Cr.nc.z_6.ld1.psl.v1.0.0': (150.0, 600.0),
 'Cs.nc.z_1.ld1.psl.v1.0.0': (30.0, 60.0),
 'F.nc.z_7.ld1.psl.v1.0.0': (200.0, 400.0),
 'Fe.nc.z_8.ld1.psl.v1.0.0': (200.0, 700.0),
 'Fr.nc.z_1.ld1.psl.v1.0.0': (30.0, 60.0),
 'Ga.nc.z_3.ld1.psl.v1.0.0': (45.0, 180.0),
 'H.nc.z_1.ld1.psl.v1.0.0': (150.0, 300.0),
 'He.nc.z_2.ld1.psl.v1.0.0': (200.0, 400.0),
 'Hf.nc.z_4.ld1.psl.v1.0.0': (200.0, 800.0),
 'Hg.nc.z_12.ld1.psl.v1.0.0': (200.0, 800.0),
 'I.nc.z_7.ld1.psl.v1.0.0': (200.0, 800.0),
 'In.nc.z_3.ld1.psl.v1.0.0': (45.0, 112.0),
 'Ir.nc.z_9.ld1.psl.v1.0.0': (100.0, 400.0),
 'K.nc.z_1.ld1.psl.v1.0.0': (30.0, 60.0),
 'Kr.nc.z_8.ld1.psl.v1.0.0': (200.0, 800.0),
 'Li.nc.z_1.ld1.psl.v1.0.0': (30.0, 105.0),
 'Mg.nc.z_2.ld1.psl.v1.0.0': (40.0, 160.0),
 'Mo.nc.z_6.ld1.psl.v1.0.0': (65.0, 260.0),
 'N.nc.z_5.ld1.psl.v1.0.0': (150.0, 300.0),
 'Na.nc.z_1.ld1.psl.v1.0.0': (30.0, 60.0),
 'Nb.nc.z_5.ld1.psl.v1.0.0': (40.0, 140.0),
 'Ne.nc.z_8.ld1.psl.v1.0.0': (200.0, 800.0),
 'Ni.nc.z_10.ld1.psl.v1.0.0': (200.0, 800.0),
 'O.nc.z_6.ld1.psl.v1.0.0': (150.0, 420.0),
 'Os.nc.z_8.ld1.psl.v1.0.0': (150.0, 525.0),
 'P.nc.z_5.ld1.psl.v1.0.0': (150.0, 450.0),
 'Pb.nc.z_4.ld1.psl.v1.0.0': (75.0, 300.0),
 'Pd.nc.z_10.ld1.psl.v1.0.0': (150.0, 525.0),
 'Pt.nc.z_10.ld1.psl.v1.0.0': (200.0, 700.0),
 'Ra.nc.z_2.ld1.psl.v1.0.0': (80.0, 200.0),
 'Rb.nc.z_1.ld1.psl.v1.0.0': (30.0, 60.0),
 'Re.nc.z_7.ld1.psl.v1.0.0': (120.0, 360.0),
 'Rh.nc.z_9.ld1.psl.v1.0.0': (100.0, 400.0),
 'Rn.nc.z_8.ld1.psl.v1.0.0': (200.0, 800.0),
 'Ru.nc.z_8.ld1.psl.v1.0.0': (75.0, 225.0),
 'S.nc.z_6.ld1.psl.v1.0.0': (150.0, 600.0),
 'Sb.nc.z_5.ld1.psl.v1.0.0': (100.0, 350.0),
 'Sc.nc.z_3.ld1.psl.v1.0.0': (40.0, 140.0),
 'Se.nc.z_6.ld1.psl.v1.0.0': (200.0, 800.0),
 'Si.nc.z_4.ld1.psl.v1.0.0': (80.0, 320.0),
 'Sn.nc.z_4.ld1.psl.v1.0.0': (90.0, 360.0),
 'Sr.nc.z_2.ld1.psl.v1.0.0': (30.0, 60.0),
 'Tc.nc.z_7.ld1.psl.v1.0.0': (70.0, 175.0),
 'Te.nc.z_6.ld1.psl.v1.0.0': (200.0, 700.0),
 'Ti.nc.z_4.ld1.psl.v1.0.0': (75.0, 300.0),
 'Tl.nc.z_3.ld1.psl.v1.0.0': (200.0, 800.0),
 'V.nc.z_5.ld1.psl.v1.0.0': (120.0, 260.0),
 'Xe.nc.z_8.ld1.psl.v1.0.0': (200.0, 800.0),
 'Y.nc.z_3.ld1.psl.v1.0.0': (30.0, 60.0),
 'Zn.nc.z_12.ld1.psl.v1.0.0': (150.0, 300.0),
 'Zr.nc.z_4.ld1.psl.v1.0.0': (35.0, 88.0),
}

comment = "nc-ld1-v1.0.0-dojo-oxygen"
computer = 'eiger-mc-mr33-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 16 for eiger
base_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/NC-PSL1.0.0'
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
