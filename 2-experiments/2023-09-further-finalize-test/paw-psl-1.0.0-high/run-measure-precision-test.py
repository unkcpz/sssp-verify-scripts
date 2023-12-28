import os

pseudos = {
 'Ag.paw.z_19.ld1.psl.v1.0.0-high': (100.0, 675.0),
 'Al.paw.z_3.ld1.psl.v1.0.0-high': (30.0, 180.0),
 'Ar.paw.z_8.ld1.psl.v1.0.0-high': (120.0, 720.0),
 'As.paw.z_15.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'At.paw.z_17.ld1.psl.v1.0.0-high': (120.0, 750.0),
 'Au.paw.z_33.ld1.psl.v1.0.0-high': (100.0, 600.0),
 'B.paw.z_3.ld1.psl.v1.0.0-high': (55.0, 440.0),
 'Ba.paw.z_10.ld1.psl.v1.0.0-high': (40.0, 240.0),
 'Be.paw.z_4.ld1.psl.v1.0.0-high': (120.0, 960.0),
 'Bi.paw.z_15.ld1.psl.v1.0.0-high': (75.0, 560.0),
 'Br.paw.z_17.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'C.paw.z_4.ld1.psl.v1.0.0-high': (60.0, 360.0),
 'Ca.paw.z_10.ld1.psl.v1.0.0-high': (60.0, 360.0),
 'Cd.paw.z_20.ld1.psl.v1.0.0-high': (120.0, 840.0),
 'Cl.paw.z_7.ld1.psl.v1.0.0-high': (35.0, 210.0),
 'Cr.paw.z_14.ld1.psl.v1.0.0-high': (150.0, 2700.0),
 'Cs.paw.z_9.ld1.psl.v1.0.0-high': (35.0, 210.0),
 'Cu.paw.z_19.ld1.psl.v1.0.0-high': (120.0, 720.0),
 'F.paw.z_7.ld1.psl.v1.0.0-high': (75.0, 450.0),
 'Fr.paw.z_19.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 'Ga.paw.z_13.ld1.psl.v1.0.0-high': (75.0, 450.0),
 'Ge.paw.z_14.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'H.paw.z_1.ld1.psl.v1.0.0-high': (45.0, 270.0),
 'He.paw.z_2.ld1.psl.v1.0.0-high': (150.0, 900.0),
 'Hf.paw.z_26.ld1.psl.v1.0.0-high.spfn': (90.0, 1620.0),
 'Hf.paw.z_36.ld1.psl.v1.0.0-high.spdfn': (200.0, 1600.0),
 'Hg.paw.z_20.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 'I.paw.z_17.ld1.psl.v1.0.0-high': (100.0, 600.0),
 'In.paw.z_13.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'Ir.paw.z_31.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'K.paw.z_9.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'Kr.paw.z_18.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 'Li.paw.z_3.ld1.psl.v1.0.0-high': (120.0, 720.0),
 'Mg.paw.z_10.ld1.psl.v1.0.0-high': (120.0, 750.0),
 'Mo.paw.z_14.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'N.paw.z_5.ld1.psl.v1.0.0-high': (45.0, 270.0),
 'Na.paw.z_9.ld1.psl.v1.0.0-high': (150.0, 900.0),
 'Nb.paw.z_13.ld1.psl.v1.0.0-high': (70.0, 420.0),
 'Ne.paw.z_8.ld1.psl.v1.0.0-high': (200.0, 1200.0),
 'Ni.paw.z_18.ld1.psl.v1.0.0-high': (100.0, 800.0),
 'O.paw.z_6.ld1.psl.v1.0.0-high': (75.0, 600.0),
 'Os.paw.z_30.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'P.paw.z_5.ld1.psl.v1.0.0-high': (35.0, 210.0),
 'Pb.paw.z_14.ld1.psl.v1.0.0-high': (75.0, 450.0),
 'Pd.paw.z_18.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'Po.paw.z_16.ld1.psl.v1.0.0-high': (100.0, 750.0),
 'Pt.paw.z_32.ld1.psl.v1.0.0-high': (100.0, 600.0),
 'Ra.paw.z_20.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 'Rb.paw.z_9.ld1.psl.v1.0.0-high': (70.0, 420.0),
 'Re.paw.z_29.ld1.psl.v1.0.0-high': (70.0, 420.0),
 'Rh.paw.z_17.ld1.psl.v1.0.0-high': (90.0, 675.0),
 'Rn.paw.z_18.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 'Ru.paw.z_16.ld1.psl.v1.0.0-high': (65.0, 480.0),
 'S.paw.z_6.ld1.psl.v1.0.0-high': (35.0, 210.0),
 'Sb.paw.z_15.ld1.psl.v1.0.0-high': (75.0, 450.0),
 'Sc.paw.z_11.ld1.psl.v1.0.0-high': (80.0, 600.0),
 'Se.paw.z_16.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'Si.paw.z_4.ld1.psl.v1.0.0-high': (30.0, 180.0),
 'Sn.paw.z_14.ld1.psl.v1.0.0-high': (75.0, 450.0),
 'Sr.paw.z_10.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'Ta.paw.z_27.ld1.psl.v1.0.0-high': (70.0, 420.0),
 'Tc.paw.z_15.ld1.psl.v1.0.0-high': (60.0, 360.0),
 'Te.paw.z_16.ld1.psl.v1.0.0-high': (90.0, 630.0),
 'Ti.paw.z_12.ld1.psl.v1.0.0-high': (80.0, 525.0),
 'Tl.paw.z_13.ld1.psl.v1.0.0-high': (120.0, 720.0),
 'V.paw.z_13.ld1.psl.v1.0.0-high': (80.0, 490.0),
 'W.paw.z_28.ld1.psl.v1.0.0-high': (150.0, 1050.0),
 'Xe.paw.z_18.ld1.psl.v1.0.0-high': (200.0, 1500.0),
 'Y.paw.z_11.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'Zn.paw.z_20.ld1.psl.v1.0.0-high': (200.0, 1200.0),
 'Zr.paw.z_12.ld1.psl.v1.0.0-high': (50.0, 400.0),
}

comment = "PAW-PSL-1.0.0-high"
computer = 'eiger-mc-mr33-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 16 for eiger
base_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/PAW-PSL1.0.0-high'
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
