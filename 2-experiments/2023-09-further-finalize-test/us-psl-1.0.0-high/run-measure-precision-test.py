import os

pseudos = {'Ag.us.z_19.ld1.psl.v1.0.0-high': (100.0, 675.0),
 'Al.us.z_3.ld1.psl.v1.0.0-high': (30.0, 180.0),
 'Ar.us.z_8.ld1.psl.v1.0.0-high': (120.0, 720.0),  
 'As.us.z_15.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'At.us.z_17.ld1.psl.v1.0.0-high': (120.0, 750.0),
 'Au.us.z_33.ld1.psl.v1.0.0-high': (100.0, 600.0),
 'B.us.z_3.ld1.psl.v1.0.0-high': (55.0, 385.0),
 'Ba.us.z_10.ld1.psl.v1.0.0-high': (35.0, 210.0),
 'Be.us.z_4.ld1.psl.v1.0.0-high': (200.0, 1200.0),
 'Bi.us.z_15.ld1.psl.v1.0.0-high': (75.0, 560.0),
 'Br.us.z_17.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'C.us.z_4.ld1.psl.v1.0.0-high': (65.0, 390.0),  
 'Ca.us.z_10.ld1.psl.v1.0.0-high': (60.0, 360.0),
 'Cd.us.z_20.ld1.psl.v1.0.0-high': (120.0, 840.0),
 'Cl.us.z_7.ld1.psl.v1.0.0-high': (35.0, 210.0), 
 'Cs.us.z_9.ld1.psl.v1.0.0-high': (35.0, 210.0),  
 'Cu.us.z_19.ld1.psl.v1.0.0-high': (120.0, 720.0),
 'F.us.z_7.ld1.psl.v1.0.0-high': (75.0, 450.0),    
 'Fe.us.z_16.ld1.psl.v1.0.0-high': (90.0, 810.0),
 'Fr.us.z_19.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 'Ga.us.z_13.ld1.psl.v1.0.0-high': (75.0, 450.0),
 'Ge.us.z_14.ld1.psl.v1.0.0-high': (65.0, 390.0),  
 'H.us.z_1.ld1.psl.v1.0.0-high': (45.0, 270.0),  
 'He.us.z_2.ld1.psl.v1.0.0-high': (100.0, 600.0),
 'Hf.us.z_26.ld1.psl.v1.0.0-high.spfn': (90.0, 1620.0),
 'Hf.us.z_36.ld1.psl.v1.0.0-high.spdfn': (200.0, 1600.0),
 'Hg.us.z_20.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 'I.us.z_17.ld1.psl.v1.0.0-high': (100.0, 600.0),
 'In.us.z_13.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'Ir.us.z_31.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'K.us.z_9.ld1.psl.v1.0.0-high': (90.0, 540.0),  
 'Kr.us.z_18.ld1.psl.v1.0.0-high': (200.0, 1500.0),
 'Li.us.z_3.ld1.psl.v1.0.0-high': (120.0, 720.0), 
 'Mg.us.z_10.ld1.psl.v1.0.0-high': (120.0, 720.0),
 'Mo.us.z_14.ld1.psl.v1.0.0-high': (65.0, 390.0), 
 'N.us.z_5.ld1.psl.v1.0.0-high': (45.0, 270.0),  
 'Na.us.z_9.ld1.psl.v1.0.0-high': (150.0, 900.0),
 'Nb.us.z_13.ld1.psl.v1.0.0-high': (70.0, 420.0),  
 'Ne.us.z_8.ld1.psl.v1.0.0-high': (55.0, 330.0),
 'O.us.z_6.ld1.psl.v1.0.0-high': (75.0, 600.0),   
 'Os.us.z_30.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'P.us.z_5.ld1.psl.v1.0.0-high': (30.0, 180.0),                     
 'Pb.us.z_14.ld1.psl.v1.0.0-high': (75.0, 450.0),
'Pd.us.z_18.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'Po.us.z_16.ld1.psl.v1.0.0-high': (100.0, 750.0),
 'Pt.us.z_32.ld1.psl.v1.0.0-high': (100.0, 600.0),
 'Ra.us.z_20.ld1.psl.v1.0.0-high': (200.0, 1600.0),
 'Rb.us.z_9.ld1.psl.v1.0.0-high': (70.0, 420.0),
 'Re.us.z_29.ld1.psl.v1.0.0-high': (70.0, 420.0),
 'Rh.us.z_17.ld1.psl.v1.0.0-high': (90.0, 675.0),
 'Rn.us.z_18.ld1.psl.v1.0.0-high': (200.0, 1500.0),
 'Ru.us.z_16.ld1.psl.v1.0.0-high': (65.0, 480.0),
 'S.us.z_6.ld1.psl.v1.0.0-high': (35.0, 210.0),
 'Sb.us.z_15.ld1.psl.v1.0.0-high': (75.0, 450.0),
 'Sc.us.z_11.ld1.psl.v1.0.0-high': (80.0, 600.0),
 'Se.us.z_16.ld1.psl.v1.0.0-high': (90.0, 540.0),
 'Si.us.z_4.ld1.psl.v1.0.0-high': (30.0, 180.0),
 'Sn.us.z_14.ld1.psl.v1.0.0-high': (75.0, 450.0),
 'Sr.us.z_10.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'Ta.us.z_27.ld1.psl.v1.0.0-high': (70.0, 450.0),
 'Tc.us.z_15.ld1.psl.v1.0.0-high': (60.0, 400.0),
 'Te.us.z_16.ld1.psl.v1.0.0-high': (100.0, 600.0),
 'Ti.us.z_12.ld1.psl.v1.0.0-high': (80.0, 525.0),
 'Tl.us.z_13.ld1.psl.v1.0.0-high': (120.0, 720.0),
 'V.us.z_13.ld1.psl.v1.0.0-high': (120.0, 960.0),
 'W.us.z_28.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'Xe.us.z_18.ld1.psl.v1.0.0-high': (200.0, 1500.0),
 'Y.us.z_11.ld1.psl.v1.0.0-high': (65.0, 390.0),
 'Zn.us.z_20.ld1.psl.v1.0.0-high': (150.0, 900.0),
 'Zr.us.z_12.ld1.psl.v1.0.0-high': (50.0, 400.0)}

comment = "US-PSL-1.0.0-high"
computer = 'eiger-mc-mr32-mem'
#computer = 'daint-mc-mrcloud-mem'
#mpiprocs = 36 # 36 for daint
mpiprocs = 128 # 128 for eiger
#npool = 4 # 4 for daint
npool = 16 # 16 for eiger
base_path = '/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe/US-PSL1.0.0-high'
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
