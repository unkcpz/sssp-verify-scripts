import os

pseudos = {
 'Ag.nc.z_19.oncvpsp3.sg15.v1.0': (50.0, 140.0),
 'Al.paw.z_3.ld1.psl.v1.0.0-high': (30.0, 180.0),
 #'Ar.paw.z_8.ld1.psl.v1.0.0-high': (120.0, 720.0),
 #'As.nc.z_15.oncvpsp3.dojo.v0.4.1-std': (80.0, 800.0),
 #'Au.nc.z_19.oncvpsp3.sg15.v1.0': (50.0, 140.0),
 #'B.us.z_3.uspp.gbrv.v1.4': (50.0, 300.0),
 #'Ba.nc.z_10.oncvpsp4.dojo.v0.5.0': (55.0, 110.0),
 #'Be.nc.z_4.oncvpsp3.sg15.v1.0': (60.0, 140.0),
 #'Bi.us.z_15.uspp.gbrv.v1': (45.0, 270.0),
 #'Br.us.z_7.uspp.gbrv.v1.4': (30.0, 180.0),
 #'C.paw.z_4.ld1.psl.v1.0.0-high': (60.0, 360.0),
 #'Ca.us.z_10.uspp.gbrv.v1': (30.0, 180.0),
 #'Cd.paw.z_20.ld1.psl.v1.0.0-high': (120.0, 840.0),
 #'Cl.us.z_7.ld1.psl.v1.0.0-high': (35.0, 210.0),
 #'Co.us.z_17.uspp.gbrv.v1.2': (150.0, 2400.0),
 #'Cr.us.z_14.uspp.gbrv.v1.5': (90.0, 1440.0),
 #'Cs.nc.z_9.oncvpsp3.dojo.v0.4.1-str': (60.0, 120.0),
 #'Cu.paw.z_11.ld1.psl.v1.0.0-low': (70.0, 455.0),
 #'F.nc.z_7.oncvpsp3.dojo.v0.4.1-std': (80.0, 260.0),
 #'Fe.paw.z_16.ld1.psl.v0.2.1': (120.0, 960.0),
 #'Ga.paw.z_13.ld1.psl.v1.0.0-high': (75.0, 450.0),
 #'Ge.us.z_14.uspp.gbrv.v1.4': (40.0, 240.0),
 #'H.nc.z_1.oncvpsp3.sg15.v1.0': (70.0, 140.0),
 #'He.nc.z_2.oncvpsp3.sg15.v1.0': (80.0, 320.0),
 #'Hf.nc.z_12.oncvpsp3.dojo.v0.4.1-std': (55.0, 110.0),
 #'Hg.us.z_12.uspp.gbrv.v1': (150.0, 900.0),
 #'I.nc.z_17.oncvpsp4.sg15.v0': (65.0, 160.0),
 #'In.us.z_13.ld1.psl.v0.2.2': (50.0, 300.0),
 #'Ir.us.z_31.ld1.psl.v1.0.0-high': (90.0, 540.0),
 #'K.paw.z_9.ld1.psl.v1.0.0-high': (90.0, 540.0),
 #'Kr.paw.z_18.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 #'Li.us.z_3.uspp.gbrv.v1.4': (40.0, 240.0),
 #'Mg.us.z_10.uspp.gbrv.v1.4': (40.0, 240.0),
 #'Mn.us.z_15.uspp.gbrv.v1.5': (150.0, 2400.0),
 #'Mo.nc.z_14.oncvpsp3.sg15.v1.0': (55.0, 140.0),
 #'N.nc.z_5.oncvpsp3.dojo.v0.4.1-std': (75.0, 188.0),
 #'Na.paw.z_9.ld1.psl.v1.0.0-low': (150.0, 900.0),
 #'Nb.paw.z_13.ld1.psl.v0.3.0': (50.0, 300.0),
 #'Ne.paw.z_8.ld1.psl.v1.0.0-high': (200.0, 1200.0),
 #'Ni.us.z_18.uspp.gbrv.v1.4': (100.0, 1800.0),
 #'O.paw.z_6.ld1.psl.v0.1': (70.0, 560.0),
 #'Os.us.z_16.uspp.gbrv.v1.2': (40.0, 240.0),
 #'P.us.z_5.ld1.psl.v1.0.0-high': (30.0, 180.0),
 #'Pb.paw.z_14.ld1.psl.v0.2.2': (50.0, 300.0),
 #'Pd.nc.z_18.oncvpsp3.sg15.v1.0': (50.0, 140.0),
 #'Po.us.z_16.ld1.psl.v1.0.0-high': (100.0, 750.0),
 #'Pt.us.z_32.ld1.psl.v1.0.0-high': (100.0, 600.0),
 #'Rb.nc.z_9.oncvpsp3.sg15.v1.0': (30.0, 60.0),
 #'Re.us.z_15.uspp.gbrv.v1.2': (65.0, 390.0),
 #'Rh.nc.z_17.oncvpsp3.sg15.v1.0': (65.0, 130.0),
 #'Rn.paw.z_18.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 #'Ru.nc.z_16.oncvpsp3.sg15.v1.0': (60.0, 120.0),
 #'S.us.z_6.uspp.gbrv.v1.4': (30.0, 180.0),
 #'Sb.us.z_15.uspp.gbrv.v1.4': (40.0, 300.0),
 #'Sc.paw.z_11.ld1.psl.v0.2.3': (75.0, 450.0),
 #'Se.us.z_6.uspp.gbrv.v1': (35.0, 262.0),
 #'Si.us.z_4.ld1.psl.v1.0.0-high': (30.0, 180.0),
 #'Sn.us.z_14.uspp.gbrv.v1.4': (90.0, 540.0),
 #'Sr.us.z_10.uspp.gbrv.v1': (40.0, 240.0),
 #'Ta.us.z_13.uspp.gbrv.v1': (50.0, 300.0),
 #'Tc.nc.z_15.oncvpsp3.sg15.v1.0': (60.0, 120.0),
 #'Te.us.z_6.ld1.psl.v1.0.0-low': (30.0, 180.0),
 #'Ti.us.z_12.uspp.gbrv.v1.4': (35.0, 210.0),
 #'Tl.us.z_13.uspp.gbrv.v1.2': (90.0, 540.0),
 #'V.us.z_13.uspp.gbrv.v1.4': (75.0, 525.0),
 #'W.us.z_14.uspp.gbrv.v1.2': (50.0, 300.0),
 #'Xe.paw.z_18.ld1.psl.v1.0.0-high': (200.0, 1500.0),
 #'Y.us.z_11.uspp.gbrv.v1.4': (45.0, 270.0),
 #'Zn.us.z_20.uspp.gbrv.v1': (200.0, 1200.0),
 #'Zr.us.z_12.uspp.gbrv.v1': (30.0, 180.0),
}

comment = "SSSP-1.3-precision"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud-mem'
mpiprocs = 36 # 128 for eiger
npool = 4 # 8 for eiger
base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/MIX-SSSP-precision-1.3.0-recollected"

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = os.path.join(base_path, pseudo)
    command = f"aiida-sssp-workflow launch --property measure.precision --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    # print(command)
