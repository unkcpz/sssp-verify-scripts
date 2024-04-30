from pathlib import Path

pseudos = {
 'Ag.us.z_19.uspp.gbrv.v1.4': (120.0, 720.0),
 'Al.us.z_3.uspp.gbrv.v1': (30.0, 180.0),
 'As.us.z_5.uspp.gbrv.v1': (60.0, 360.0),
 'Au.us.z_11.uspp.gbrv.v1': (65.0, 488.0),
 'B.us.z_3.uspp.gbrv.v1.4': (50.0, 300.0),
 'Ba.us.z_10.uspp.gbrv.v1': (45.0, 270.0),
 'Be.us.z_4.uspp.gbrv.v1.4': (45.0, 270.0),
 'Bi.us.z_15.uspp.gbrv.v1': (45.0, 270.0),
 'Br.us.z_7.uspp.gbrv.v1.4': (30.0, 180.0),
 'C.us.z_4.uspp.gbrv.v1.2': (200.0, 1200.0),
 'Ca.us.z_10.uspp.gbrv.v1': (30.0, 180.0),
 'Cd.us.z_12.uspp.gbrv.v1': (150.0, 900.0),
 'Cl.us.z_7.uspp.gbrv.v1.4': (35.0, 210.0),
 'Co.us.z_17.uspp.gbrv.v1.2': (150.0, 2400.0),
 'Cr.us.z_14.uspp.gbrv.v1.5': (200.0, 3600.0),
 'Cs.us.z_9.uspp.gbrv.v1': (35.0, 210.0),
 'Cu.us.z_19.uspp.gbrv.v1.2': (200.0, 1200.0),
 'F.us.z_7.uspp.gbrv.v1.4': (45.0, 270.0),
 'Fe.us.z_16.uspp.gbrv.v1.5': (150.0, 2400.0),
 'Ga.us.z_19.uspp.gbrv.v1.4': (200.0, 1200.0),
 'Ge.us.z_14.uspp.gbrv.v1.4': (40.0, 240.0),
 'H.us.z_1.uspp.gbrv.v1.4': (30.0, 180.0),
 'Hf.us.z_12.uspp.gbrv.plus4_v1': (200.0, 1600.0),
 'Hf.us.z_12.uspp.gbrv.v1': (200.0, 3200.0),
 'Hg.us.z_12.uspp.gbrv.v1': (150.0, 1200.0),
 'I.us.z_7.uspp.gbrv.v1': (30.0, 210.0),
 'In.us.z_13.uspp.gbrv.v1.4': (150.0, 900.0),
 'Ir.us.z_15.uspp.gbrv.v1.2': (90.0, 720.0),
 'K.us.z_9.uspp.gbrv.v1.4': (80.0, 560.0),
 'Li.us.z_3.uspp.gbrv.v1.4': (40.0, 240.0),
 'Mg.us.z_10.uspp.gbrv.v1.4': (40.0, 240.0),
 'Mn.us.z_15.uspp.gbrv.v1.5': (150.0, 2400.0),
 'Mo.us.z_14.uspp.gbrv.v1': (50.0, 375.0),
 'N.us.z_5.uspp.gbrv.v1.2': (55.0, 330.0),
 'Na.us.z_9.uspp.gbrv.v1.5': (200.0, 1200.0),
 'Nb.us.z_13.uspp.gbrv.v1': (50.0, 375.0),
 'Ni.us.z_18.uspp.gbrv.v1.4': (100.0, 1800.0),
 'O.us.z_6.uspp.gbrv.v1.2': (55.0, 990.0),
 'Os.us.z_16.uspp.gbrv.v1.2': (40.0, 240.0),
 'P.us.z_5.uspp.gbrv.v1.5': (30.0, 180.0),
 'Pb.us.z_14.uspp.gbrv.v1': (35.0, 280.0),
 'Pd.us.z_16.uspp.gbrv.v1.4': (40.0, 240.0),
 'Pt.us.z_16.uspp.gbrv.v1.4': (65.0, 488.0),
 'Rb.us.z_9.uspp.gbrv.v1': (30.0, 180.0),
 'Re.us.z_15.uspp.gbrv.v1.2': (65.0, 390.0),
 'Rh.us.z_15.uspp.gbrv.v1.4': (45.0, 270.0),
 'Ru.us.z_16.uspp.gbrv.v1.2': (50.0, 300.0),
 'S.us.z_6.uspp.gbrv.v1.4': (30.0, 180.0),
 'Sb.us.z_15.uspp.gbrv.v1.4': (40.0, 300.0),
 'Sc.us.z_11.uspp.gbrv.v1': (45.0, 270.0),
 'Se.us.z_6.uspp.gbrv.v1': (35.0, 262.0),
 'Si.us.z_4.uspp.gbrv.v1': (30.0, 180.0),
 'Sn.us.z_14.uspp.gbrv.v1.4': (90.0, 540.0),
 'Sr.us.z_10.uspp.gbrv.v1': (40.0, 240.0),
 'Ta.us.z_13.uspp.gbrv.v1': (50.0, 300.0),
 'Tc.us.z_15.uspp.gbrv.v1': (150.0, 975.0),
 'Te.us.z_6.uspp.gbrv.v1': (200.0, 1500.0),
 'Ti.us.z_12.uspp.gbrv.v1.4': (35.0, 210.0),
 'Tl.us.z_13.uspp.gbrv.v1.2': (70.0, 520.0),
 'V.us.z_13.uspp.gbrv.v1.4': (120.0, 720.0),
 'W.us.z_14.uspp.gbrv.v1.2': (50.0, 300.0),
 'Y.us.z_11.uspp.gbrv.v1.4': (45.0, 270.0),
 'Zn.us.z_20.uspp.gbrv.v1': (200.0, 1200.0),
 'Zr.us.z_12.uspp.gbrv.v1': (30.0, 180.0),
}


comment = "US-GBRV-v1.x-UPF2"
computer = 'eiger-mem-hq'
mpiprocs = 32
npool = 4

# The base path for this repo contains all the source pesudos and scripts to run the verification
base_path = Path("/home/jyu/Project/sssp-project/sssp-verify-scripts")
lib_base_path = base_path / "libraries-pbe/US-GBRV-1.x-upf2"

sssp_oxygen = base_path / "common/pseudos/O.paw.pbe.z_6.ld1.psl.v0.1.upf"
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
