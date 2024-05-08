import os
from pathlib import Path

pseudos = {
 #'Ag.us.pbe.z_19.uspp.gbrv.v1.4': (120.0, 720.0),
 'Al.us.pbe.z_3.uspp.gbrv.v1': (30.0, 180.0),
 #'As.us.pbe.z_5.uspp.gbrv.v1': (60.0, 360.0),
 #'Au.us.pbe.z_11.uspp.gbrv.v1': (65.0, 488.0),
 #'B.us.pbe.z_3.uspp.gbrv.v1.4': (50.0, 300.0),
 #'Ba.us.pbe.z_10.uspp.gbrv.v1': (45.0, 270.0),
 #'Be.us.pbe.z_4.uspp.gbrv.v1.4': (45.0, 270.0),
 #'Bi.us.pbe.z_15.uspp.gbrv.v1': (45.0, 270.0),
 #'Br.us.pbe.z_7.uspp.gbrv.v1.4': (30.0, 180.0),
 #'C.us.pbe.z_4.uspp.gbrv.v1.2': (200.0, 1200.0),
 #'Ca.us.pbe.z_10.uspp.gbrv.v1': (30.0, 180.0),
 #'Cd.us.pbe.z_12.uspp.gbrv.v1': (150.0, 900.0),
 #'Cl.us.pbe.z_7.uspp.gbrv.v1.4': (35.0, 210.0),
 #'Co.us.pbe.z_17.uspp.gbrv.v1.2': (150.0, 2400.0),
 #'Cr.us.pbe.z_14.uspp.gbrv.v1.5': (200.0, 3600.0),
 #'Cs.us.pbe.z_9.uspp.gbrv.v1': (35.0, 210.0),
 #'Cu.us.pbe.z_19.uspp.gbrv.v1.2': (200.0, 1200.0),
 #'F.us.pbe.z_7.uspp.gbrv.v1.4': (45.0, 270.0),
 #'Fe.us.pbe.z_16.uspp.gbrv.v1.5': (150.0, 2400.0),
 #'Ga.us.pbe.z_19.uspp.gbrv.v1.4': (200.0, 1200.0),
 #'Ge.us.pbe.z_14.uspp.gbrv.v1.4': (40.0, 240.0),
 #'H.us.pbe.z_1.uspp.gbrv.v1.4': (30.0, 180.0),
 #'Hf.us.pbe.z_12.uspp.gbrv.plus4_v1': (200.0, 1600.0),
 #'Hf.us.pbe.z_12.uspp.gbrv.v1': (200.0, 3200.0),
 #'Hg.us.pbe.z_12.uspp.gbrv.v1': (150.0, 1200.0),
 #'I.us.pbe.z_7.uspp.gbrv.v1': (30.0, 210.0),
 #'In.us.pbe.z_13.uspp.gbrv.v1.4': (150.0, 900.0),
 #'Ir.us.pbe.z_15.uspp.gbrv.v1.2': (90.0, 720.0),
 #'K.us.pbe.z_9.uspp.gbrv.v1.4': (80.0, 560.0),
 #'Li.us.pbe.z_3.uspp.gbrv.v1.4': (40.0, 240.0),
 #'Mg.us.pbe.z_10.uspp.gbrv.v1.4': (40.0, 240.0),
 #'Mn.us.pbe.z_15.uspp.gbrv.v1.5': (150.0, 2400.0),
 #'Mo.us.pbe.z_14.uspp.gbrv.v1': (50.0, 375.0),
 #'N.us.pbe.z_5.uspp.gbrv.v1.2': (55.0, 330.0),
 #'Na.us.pbe.z_9.uspp.gbrv.v1.5': (200.0, 1200.0),
 #'Nb.us.pbe.z_13.uspp.gbrv.v1': (50.0, 375.0),
 #'Ni.us.pbe.z_18.uspp.gbrv.v1.4': (100.0, 1800.0),
 #'O.us.pbe.z_6.uspp.gbrv.v1.2': (55.0, 990.0),
 #'Os.us.pbe.z_16.uspp.gbrv.v1.2': (40.0, 240.0),
 #'P.us.pbe.z_5.uspp.gbrv.v1.5': (30.0, 180.0),
 #'Pb.us.pbe.z_14.uspp.gbrv.v1': (35.0, 280.0),
 #'Pd.us.pbe.z_16.uspp.gbrv.v1.4': (40.0, 240.0),
 #'Pt.us.pbe.z_16.uspp.gbrv.v1.4': (65.0, 488.0),
 #'Rb.us.pbe.z_9.uspp.gbrv.v1': (30.0, 180.0),
 #'Re.us.pbe.z_15.uspp.gbrv.v1.2': (65.0, 390.0),
 #'Rh.us.pbe.z_15.uspp.gbrv.v1.4': (45.0, 270.0),
 #'Ru.us.pbe.z_16.uspp.gbrv.v1.2': (50.0, 300.0),
 #'S.us.pbe.z_6.uspp.gbrv.v1.4': (30.0, 180.0),
 #'Sb.us.pbe.z_15.uspp.gbrv.v1.4': (40.0, 300.0),
 #'Sc.us.pbe.z_11.uspp.gbrv.v1': (45.0, 270.0),
 #'Se.us.pbe.z_6.uspp.gbrv.v1': (35.0, 262.0),
 #'Si.us.pbe.z_4.uspp.gbrv.v1': (30.0, 180.0),
 #'Sn.us.pbe.z_14.uspp.gbrv.v1.4': (90.0, 540.0),
 #'Sr.us.pbe.z_10.uspp.gbrv.v1': (40.0, 240.0),
 #'Ta.us.pbe.z_13.uspp.gbrv.v1': (50.0, 300.0),
 #'Tc.us.pbe.z_15.uspp.gbrv.v1': (150.0, 975.0),
 #'Te.us.pbe.z_6.uspp.gbrv.v1': (200.0, 1500.0),
 #'Ti.us.pbe.z_12.uspp.gbrv.v1.4': (35.0, 210.0),
 #'Tl.us.pbe.z_13.uspp.gbrv.v1.2': (70.0, 520.0),
 #'V.us.pbe.z_13.uspp.gbrv.v1.4': (120.0, 720.0),
 #'W.us.pbe.z_14.uspp.gbrv.v1.2': (50.0, 300.0),
 #'Y.us.pbe.z_11.uspp.gbrv.v1.4': (45.0, 270.0),
 #'Zn.us.pbe.z_20.uspp.gbrv.v1': (200.0, 1200.0),
 #'Zr.us.pbe.z_12.uspp.gbrv.v1': (30.0, 180.0),
}


comment = "US-GBRV-v1.x-UPF2"
computer = 'eiger-hq'
num_cpus = 32
memory_mb = 120000 # mb
npool = 4

# The base path for this repo contains all the source pesudos and scripts to run the verification
base_path = Path("/home/jyu/project/sssp-project/sssp-verify-scripts")
lib_base_path = base_path / "libraries-pbe/US-GBRV-1.x-upf2"

sssp_oxygen = base_path / "common/pseudos/O.paw.pbe.z_6.ld1.psl.v0.1.upf"

oxygen_ecutwfc = 70
oxygen_ecutrho = 560

for pseudo, (wfc_cutoff, rho_cutoff) in pseudos.items():
    pseudo = pseudo + '.upf'
    pseudo_path = lib_base_path / pseudo
    command = f"aiida-sssp-workflow launch --property measure.precision --oxygen-pseudo {sssp_oxygen} --oxygen-ecutwfc {oxygen_ecutwfc} --oxygen-ecutrho {oxygen_ecutrho} --ecutwfc {wfc_cutoff} --ecutrho {rho_cutoff} --pw-code pw-7.0@{computer} --protocol acwf --withmpi True --resources num_cpus {num_cpus} --resources memory_mb {memory_mb} --npool {npool} --no-daemon --comment {comment} -- {pseudo_path}"
    os.system(command)
    print(f"Launched {pseudo}")
