# verdi run run-convergence-wrt-diff-conf.py

import os

pseudos = [
"Ag.nc.z_19.oncvpsp3.sg15.v1.0.upf",
#"Al.paw.z_3.ld1.psl.v1.0.0-high.upf",
#"Am.paw.z_17.ld1.uni-marburg.v0.upf",
#"Ar.paw.z_8.ld1.psl.v1.0.0-high.upf",
#"As.nc.z_15.oncvpsp3.dojo.v0.4.1-std.upf",
#"Au.nc.z_19.oncvpsp3.sg15.v1.0.upf",
#"Ba.nc.z_10.oncvpsp4.dojo.v0.5.0.upf",
#"Be.nc.z_4.oncvpsp3.sg15.v1.0.upf",
#"Bi.us.z_15.uspp.gbrv.v1.upf",
#"Bk.paw.z_19.ld1.uni-marburg.v0.upf",
#"Br.us.z_7.uspp.gbrv.v1.4.upf",
#"B.us.z_3.uspp.gbrv.v1.4.upf",
#"Ca.us.z_10.uspp.gbrv.v1.upf",
#"Cd.paw.z_20.ld1.psl.v1.0.0-high.upf",
#"Cf.paw.z_20.ld1.uni-marburg.v0.upf",
#"Cl.us.z_7.ld1.psl.v1.0.0-high.upf",
#"Cm.paw.z_18.ld1.uni-marburg.v0.upf",
#"Co.us.z_17.uspp.gbrv.v1.2.upf",
#"C.paw.z_4.ld1.psl.v1.0.0-high.upf",
#"Cr.us.z_14.uspp.gbrv.v1.5.upf",
#"Cs.nc.z_9.oncvpsp3.dojo.v0.4.1-str.upf",
#"Cu.paw.z_11.ld1.psl.v1.0.0-low.upf",
#"Es.paw.z_21.ld1.uni-marburg.v0.upf",
#"Fe.paw.z_16.ld1.psl.v0.2.1.upf",
#"Fm.paw.z_22.ld1.uni-marburg.v0.upf",
#"F.nc.z_7.oncvpsp3.dojo.v0.4.1-std.upf",
#"Ga.paw.z_13.ld1.psl.v1.0.0-high.upf",
#"Ge.us.z_14.uspp.gbrv.v1.4.upf",
#"He.nc.z_2.oncvpsp3.sg15.v1.0.upf",
#"Hf.nc.z_12.oncvpsp3.dojo.v0.4.1-std.upf",
#"Hg.us.z_12.uspp.gbrv.v1.upf",
#"H.nc.z_1.oncvpsp3.sg15.v1.0.upf",
#"I.nc.z_17.oncvpsp4.sg15.v0.upf",
#"In.us.z_13.ld1.psl.v0.2.2.upf",
#"Ir.us.z_31.ld1.psl.v1.0.0-high.upf",
#"K.paw.z_9.ld1.psl.v1.0.0-high.upf",
#"Kr.paw.z_18.ld1.psl.v1.0.0-high.upf",
#"Li.us.z_3.uspp.gbrv.v1.4.upf",
#"Lr.paw.z_25.ld1.uni-marburg.v0.upf",
#"Md.paw.z_23.ld1.uni-marburg.v0.upf",
#"Mg.us.z_10.uspp.gbrv.v1.4.upf",
#"Mn.us.z_15.uspp.gbrv.v1.5.upf",
#"Mo.nc.z_14.oncvpsp3.sg15.v1.0.upf",
#"Na.paw.z_9.ld1.psl.v1.0.0-low.upf",
#"Nb.paw.z_13.ld1.psl.v0.3.0.upf",
#"Ne.paw.z_8.ld1.psl.v1.0.0-high.upf",
#"Ni.us.z_18.uspp.gbrv.v1.4.upf",
#"N.nc.z_5.oncvpsp3.dojo.v0.4.1-std.upf",
#"No.paw.z_24.ld1.uni-marburg.v0.upf",
#"Np.paw.z_15.ld1.uni-marburg.v0.upf",
#"O.paw.z_6.ld1.psl.v0.1.upf",
#"Os.us.z_16.uspp.gbrv.v1.2.upf",
#"Pa.paw.z_13.ld1.uni-marburg.v0.upf",
#"Pb.paw.z_14.ld1.psl.v0.2.2.upf",
#"Pd.nc.z_18.oncvpsp3.sg15.v1.0.upf",
#"Po.us.z_16.ld1.psl.v1.0.0-high.upf",
#"Pt.us.z_32.ld1.psl.v1.0.0-high.upf",
#"Pu.paw.z_16.ld1.uni-marburg.v0.upf",
#"P.us.z_5.ld1.psl.v1.0.0-high.upf",
#"Rb.nc.z_9.oncvpsp3.sg15.v1.0.upf",
#"Re.us.z_15.uspp.gbrv.v1.2.upf",
#"Rh.nc.z_17.oncvpsp3.sg15.v1.0.upf",
#"Rn.paw.z_18.ld1.psl.v1.0.0-high.upf",
#"Ru.nc.z_16.oncvpsp3.sg15.v1.0.upf",
#"Sb.us.z_15.uspp.gbrv.v1.4.upf",
#"Sc.paw.z_11.ld1.psl.v0.2.3.upf",
#"Se.us.z_6.uspp.gbrv.v1.upf",
#"Si.us.z_4.ld1.psl.v1.0.0-high.upf",
#"Sn.us.z_14.uspp.gbrv.v1.4.upf",
#"Sr.us.z_10.uspp.gbrv.v1.upf",
#"S.us.z_6.uspp.gbrv.v1.4.upf",
#"Ta.us.z_13.uspp.gbrv.v1.upf",
#"Tc.nc.z_15.oncvpsp3.sg15.v1.0.upf",
#"Te.us.z_6.ld1.psl.v1.0.0-low.upf",
#"Th.paw.z_12.ld1.uni-marburg.v0.upf",
#"Ti.us.z_12.uspp.gbrv.v1.4.upf",
#"Tl.us.z_13.uspp.gbrv.v1.2.upf",
#"U.paw.z_14.ld1.uni-marburg.v0.upf",
#"V.us.z_13.uspp.gbrv.v1.4.upf",
#"W.us.z_14.uspp.gbrv.v1.2.upf",
#"Xe.paw.z_18.ld1.psl.v1.0.0-high.upf",
#"Y.us.z_11.uspp.gbrv.v1.4.upf",
#"Zn.us.z_20.uspp.gbrv.v1.upf",
#"Zr.us.z_12.uspp.gbrv.v1.upf",
]

criteria = "precision"
#computer = 'eiger-mc-mr32-mem'
computer = 'daint-mc-mrcloud'
mpiprocs = 32 # 128 for eiger
npool = 4 # 8 for eiger
base_path = "/home/jyu/Projects/WP-SSSP/sssp-verify-scripts/libraries-pbe"
p_path = "MIX-SSSP-precision-1.3.0-recollected"
comment = "MIX-SSSP-1.3-precision-recollected"

for pseudo in pseudos:
    pseudo_path = os.path.join(base_path, p_path, pseudo)
    command = f"aiida-sssp-workflow launch --property convergence --pw-code pw-7.0@{computer} --ph-code ph-7.0@{computer} --protocol acwf --cutoff-control standard --criteria {criteria} --withmpi True --num-mpiprocs {mpiprocs} --npool {npool} --comment {comment} -- {pseudo_path}"
    os.system(command)
    # print(command)
    print(f"Launched {pseudo}")
