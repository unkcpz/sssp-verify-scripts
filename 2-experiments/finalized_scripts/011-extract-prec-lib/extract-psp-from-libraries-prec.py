import os
from pathlib import Path
import shutil

psp_pbe = [
    "Ac.us.pbe.z_11.ld1.psl.v1.0.0-high.upf",
    "Ag.nc.pbe.z_19.oncvpsp4.spms.v1.upf",
    "Al.us.pbe.z_3.ld1.psl.v1.0.0-low.upf",
    "Am.us.pbe.z_17.ld1.psl.v1.0.0-high.upf",
    "Ar.paw.pbe.z_8.ld1.psl.v1.0.0-low.upf",
    "As.paw.pbe.z_15.atompaw.jth.v1.1-std.upf",
    "At.us.pbe.z_17.ld1.psl.v1.0.0-high.upf",
    "Au.nc.pbe.z_19.oncvpsp4.spms.v1.upf",
    "B.paw.pbe.z_3.atompaw.jth.v1.1-std.upf",
    "Ba.nc.pbe.z_10.oncvpsp4.dojo.v0.5.0-std.upf",
    "Be.us.pbe.z_2.ld1.psl.v1.0.0-low.n.upf",
    "Bi.us.pbe.z_15.uspp.gbrv.v1.upf",
    "Br.nc.pbe.z_7.oncvpsp4.spms.v1.upf",
    "C.us.pbe.z_4.uspp.gbrv.v1.2.upf",
    "Ca.us.pbe.z_10.uspp.gbrv.v1.upf",
    "Cd.paw.pbe.z_12.atompaw.jth.v1.1-std.upf",
    "Ce.paw.pbe.z_12.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Cl.us.pbe.z_7.ld1.psl.v1.0.0-high.upf",
    "Co.nc.pbe.z_17.oncvpsp4.spms.v1.upf",
    "Cr.us.pbe.z_14.uspp.gbrv.v1.5.upf",
    "Cs.nc.pbe.z_9.oncvpsp4.spms.v1.upf",
    "Cu.paw.pbe.z_11.ld1.psl.v1.0.0-low.upf",
    "Dy.paw.pbe.z_20.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Er.paw.pbe.z_22.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Eu.paw.pbe.z_17.atompaw.wentzcovitch.v1.0.legacy.upf",
    "F.us.pbe.z_7.ld1.psl.v0.1.upf",
    "Fe.paw.pbe.z_16.ld1.psl.v0.2.1.upf",
    "Fr.paw.pbe.z_19.ld1.psl.v1.0.0-high.upf",
    "Ga.paw.pbe.z_13.atompaw.jth.v1.1-std.upf",
    "Gd.paw.pbe.z_18.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Ge.us.pbe.z_14.uspp.gbrv.v1.4.upf",
    "H.us.pbe.z_1.ld1.psl.v1.0.0-high.upf",
    "He.paw.pbe.z_2.ld1.psl.v1.0.0-high.upf",
    "Hf.nc.pbe.z_26.oncvpsp4.sg15.v0.upf",
    "Hg.us.pbe.z_12.uspp.gbrv.v1.upf",
    "Ho.paw.pbe.z_21.atompaw.wentzcovitch.v1.0.legacy.upf",
    "I.us.pbe.z_7.uspp.gbrv.v1.upf",
    "In.us.pbe.z_13.ld1.psl.v0.2.2.upf",
    "Ir.us.pbe.z_15.uspp.gbrv.v1.2.upf",
    "K.nc.pbe.z_9.oncvpsp4.sg15.v0.upf",
    "La.paw.pbe.z_11.ld1.psl.v1.0.0-high.spfn.upf",
    "Li.us.pbe.z_3.uspp.gbrv.v1.4.upf",
    "Lu.paw.pbe.z_25.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Mg.us.pbe.z_10.uspp.gbrv.v1.4.upf",
    "Mn.us.pbe.z_15.uspp.gbrv.v1.5.upf",
    "Mo.paw.pbe.z_14.atompaw.jth.v1.1-std.upf",
    "N.paw.pbe.z_5.ld1.psl.v0.1.upf",
    "Na.paw.pbe.z_9.ld1.psl.v0.2.upf",
    "Nb.paw.pbe.z_13.atompaw.jth.v1.1-std.upf",
    "Nd.paw.pbe.z_14.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Ne.paw.pbe.z_8.ld1.psl.v1.0.0-high.upf",
    "Ni.nc.pbe.z_18.oncvpsp4.spms.v1.upf",
    "Np.us.pbe.z_15.ld1.psl.v1.0.0-high.upf",
    "O.paw.pbe.z_6.ld1.psl.v0.1.upf",
    "Os.us.pbe.z_16.uspp.gbrv.v1.2.upf",
    "P.us.pbe.z_5.ld1.psl.v1.0.0-high.upf",
    "Pa.paw.pbe.z_13.ld1.uni-marburg.v0.upf",
    "Pb.nc.pbe.z_14.oncvpsp4.sg15.v0.upf",
    "Pd.nc.pbe.z_18.oncvpsp4.sg15.v0.upf",
    "Pm.paw.pbe.z_15.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Po.nc.pbe.z_16.oncvpsp3.dojo.v0.4.1-std.upf",
    "Pr.paw.pbe.z_13.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Pt.paw.pbe.z_18.ld1.psl.v1.0.0-low.spn.upf",
    "Pu.us.pbe.z_16.ld1.psl.v1.0.0-high.upf",
    "Ra.paw.pbe.z_20.ld1.psl.v1.0.0-high.upf",
    "Rb.us.pbe.z_9.uspp.gbrv.v1.upf",
    "Re.us.pbe.z_15.uspp.gbrv.v1.2.upf",
    "Rh.nc.pbe.z_17.oncvpsp4.spms.v1.upf",
    "Rn.paw.pbe.z_18.ld1.psl.v1.0.0-high.upf",
    "Ru.nc.pbe.z_16.oncvpsp4.sg15.v0.upf",
    "S.nc.pbe.z_6.oncvpsp4.spms.v1.upf",
    "Sb.nc.pbe.z_15.oncvpsp4.spms.v1.upf",
    "Sc.us.pbe.z_11.uspp.gbrv.v1.upf",
    "Se.us.pbe.z_6.uspp.gbrv.v1.upf",
    "Si.us.pbe.z_4.ld1.psl.v1.0.0-high.upf",
    "Sm.us.pbe.z_26.ld1.psl.v1.0.0-high.upf",
    "Sn.paw.pbe.z_14.atompaw.jth.v1.1-std.upf",
    "Sr.us.pbe.z_10.uspp.gbrv.v1.upf",
    "Ta.nc.pbe.z_13.oncvpsp4.spms.v1.upf",
    "Tb.paw.pbe.z_19.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Tc.us.pbe.z_15.uspp.gbrv.v1.upf",
    "Te.us.pbe.z_6.ld1.psl.v1.0.0-low.upf",
    "Th.us.pbe.z_12.ld1.psl.v1.0.0-high.upf",
    "Ti.us.pbe.z_12.uspp.gbrv.v1.4.upf",
    "Tl.nc.pbe.z_13.oncvpsp4.sg15.v0.upf",
    "Tm.paw.pbe.z_23.atompaw.wentzcovitch.v1.0.legacy.upf",
    "U.paw.pbe.z_14.ld1.uni-marburg.v0.upf",
    "V.us.pbe.z_13.uspp.gbrv.v1.4.upf",
    "W.nc.pbe.z_14.oncvpsp4.spms.v1.upf",
    "Xe.nc.pbe.z_8.oncvpsp4.dojo.v0.5.0-std.upf",
    "Y.us.pbe.z_11.uspp.gbrv.v1.4.upf",
    "Yb.paw.pbe.z_24.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Zn.paw.pbe.z_12.atompaw.jth.v1.1-std.upf",
    "Zr.us.pbe.z_12.uspp.gbrv.v1.upf",
]

psp_pbesol = [
    "Ac.us.pbesol.z_11.ld1.psl.v1.0.0-high.upf",
    "Ag.nc.pbesol.z_19.oncvpsp4.spms.v1.upf",
    "Al.us.pbesol.z_3.ld1.psl.v1.0.0-low.upf",
    "Am.us.pbesol.z_17.ld1.psl.v1.0.0-high.upf",
    "Ar.paw.pbesol.z_8.ld1.psl.v1.0.0-low.upf",
    "As.paw.pbesol.z_15.atompaw.jth.v1.1-std.upf",
    "At.us.pbesol.z_17.ld1.psl.v1.0.0-high.upf",
    "Au.nc.pbesol.z_19.oncvpsp4.spms.v1.upf",
    "B.paw.pbesol.z_3.atompaw.jth.v1.1-std.upf",
    "Ba.nc.pbesol.z_10.oncvpsp4.dojo.v0.5.0-std.upf",
    "Be.us.pbesol.z_2.ld1.psl.v1.0.0-low.n.upf",
    "Bi.us.pbesol.z_15.uspp.gbrv.v1.upf",
    "Br.nc.pbesol.z_7.oncvpsp4.spms.v1.upf",
    "C.us.pbesol.z_4.uspp.gbrv.v1.2.upf",
    "Ca.us.pbesol.z_10.uspp.gbrv.v1.upf",
    "Cd.paw.pbesol.z_12.atompaw.jth.v1.1-std.upf",
    "Ce.paw.pbesol.z_12.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Cl.us.pbesol.z_7.ld1.psl.v1.0.0-high.upf",
    "Co.nc.pbesol.z_17.oncvpsp4.spms.v1.upf",
    "Cr.us.pbesol.z_14.uspp.gbrv.v1.5.upf",
    "Cs.nc.pbesol.z_9.oncvpsp4.spms.v1.upf",
    "Cu.paw.pbesol.z_11.ld1.psl.v1.0.0-low.upf",
    "Dy.paw.pbesol.z_20.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Er.paw.pbesol.z_22.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Eu.paw.pbesol.z_17.atompaw.wentzcovitch.v1.0.legacy.upf",
    "F.us.pbesol.z_7.ld1.psl.v0.1.upf",
    "Fe.paw.pbesol.z_16.ld1.psl.v0.2.1.upf",
    "Fr.paw.pbesol.z_19.ld1.psl.v1.0.0-high.upf",
    "Ga.paw.pbesol.z_13.atompaw.jth.v1.1-std.upf",
    "Gd.paw.pbesol.z_18.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Ge.us.pbesol.z_14.uspp.gbrv.v1.4.upf",
    "H.us.pbesol.z_1.ld1.psl.v1.0.0-high.upf",
    "He.paw.pbesol.z_2.ld1.psl.v1.0.0-high.upf",
    "Hf.nc.pbesol.z_26.oncvpsp4.sg15.v0.upf",
    "Hg.us.pbesol.z_12.uspp.gbrv.v1.upf",
    "Ho.paw.pbesol.z_21.atompaw.wentzcovitch.v1.0.legacy.upf",
    "I.us.pbesol.z_7.uspp.gbrv.v1.upf",
    "In.us.pbesol.z_13.ld1.psl.v0.2.2.upf",
    "Ir.us.pbesol.z_15.uspp.gbrv.v1.2.upf",
    "K.nc.pbesol.z_9.oncvpsp4.sg15.v0.upf",
    "La.paw.pbesol.z_11.ld1.psl.v1.0.0-high.spfn.upf",
    "Li.us.pbesol.z_3.uspp.gbrv.v1.4.upf",
    "Lu.paw.pbesol.z_25.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Mg.us.pbesol.z_10.uspp.gbrv.v1.4.upf",
    "Mn.us.pbesol.z_15.uspp.gbrv.v1.5.upf",
    "Mo.paw.pbesol.z_14.atompaw.jth.v1.1-std.upf",
    "N.paw.pbesol.z_5.ld1.psl.v0.1.upf",
    "Na.paw.pbesol.z_9.ld1.psl.v0.2.upf",
    "Nb.paw.pbesol.z_13.atompaw.jth.v1.1-std.upf",
    "Nd.paw.pbesol.z_14.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Ne.paw.pbesol.z_8.ld1.psl.v1.0.0-high.upf",
    "Ni.nc.pbesol.z_18.oncvpsp4.spms.v1.upf",
    "Np.us.pbesol.z_15.ld1.psl.v1.0.0-high.upf",
    "O.paw.pbesol.z_6.ld1.psl.v0.1.upf",
    "Os.us.pbesol.z_16.uspp.gbrv.v1.2.upf",
    "P.us.pbesol.z_5.ld1.psl.v1.0.0-high.upf",
    "Pa.paw.pbesol.z_13.ld1.uni-marburg.v0.upf",
    "Pb.nc.pbesol.z_14.oncvpsp4.sg15.v0.upf",
    "Pd.nc.pbesol.z_18.oncvpsp4.sg15.v0.upf",
    "Pm.paw.pbesol.z_15.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Po.nc.pbesol.z_16.oncvpsp3.dojo.v0.4.1-std.upf",
    "Pr.paw.pbesol.z_13.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Pt.paw.pbesol.z_18.ld1.psl.v1.0.0-low.spn.upf",
    "Pu.us.pbesol.z_16.ld1.psl.v1.0.0-high.upf",
    "Ra.paw.pbesol.z_20.ld1.psl.v1.0.0-high.upf",
    "Rb.us.pbesol.z_9.uspp.gbrv.v1.upf",
    "Re.us.pbesol.z_15.uspp.gbrv.v1.2.upf",
    "Rh.nc.pbesol.z_17.oncvpsp4.spms.v1.upf",
    "Rn.paw.pbesol.z_18.ld1.psl.v1.0.0-high.upf",
    "Ru.nc.pbesol.z_16.oncvpsp4.sg15.v0.upf",
    "S.nc.pbesol.z_6.oncvpsp4.spms.v1.upf",
    "Sb.nc.pbesol.z_15.oncvpsp4.spms.v1.upf",
    "Sc.us.pbesol.z_11.uspp.gbrv.v1.upf",
    "Se.us.pbesol.z_6.uspp.gbrv.v1.upf",
    "Si.us.pbesol.z_4.ld1.psl.v1.0.0-high.upf",
    "Sm.us.pbesol.z_26.ld1.psl.v1.0.0-high.upf",
    "Sn.paw.pbesol.z_14.atompaw.jth.v1.1-std.upf",
    "Sr.us.pbesol.z_10.uspp.gbrv.v1.upf",
    "Ta.nc.pbesol.z_13.oncvpsp4.spms.v1.upf",
    "Tb.paw.pbesol.z_19.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Tc.us.pbesol.z_15.uspp.gbrv.v1.upf",
    "Te.us.pbesol.z_6.ld1.psl.v1.0.0-low.upf",
    "Th.us.pbesol.z_12.ld1.psl.v1.0.0-high.upf",
    "Ti.us.pbesol.z_12.uspp.gbrv.v1.4.upf",
    "Tl.nc.pbesol.z_13.oncvpsp4.sg15.v0.upf",
    "Tm.paw.pbesol.z_23.atompaw.wentzcovitch.v1.0.legacy.upf",
    "U.paw.pbesol.z_14.ld1.uni-marburg.v0.upf",
    "V.us.pbesol.z_13.uspp.gbrv.v1.4.upf",
    "W.nc.pbesol.z_14.oncvpsp4.spms.v1.upf",
    "Xe.nc.pbesol.z_8.oncvpsp4.dojo.v0.5.0-std.upf",
    "Y.us.pbesol.z_11.uspp.gbrv.v1.4.upf",
    "Yb.paw.pbesol.z_24.atompaw.wentzcovitch.v1.0.legacy.upf",
    "Zn.paw.pbesol.z_12.atompaw.jth.v1.1-std.upf",
    "Zr.us.pbesol.z_12.uspp.gbrv.v1.upf",
]


def extract(t):
    if t == "pbe":
        ll = psp_pbe
        base_folder = "../../../libraries-pbe/"
        sssp_path = "./mix-sssp-prec-pbe-lib-v2/library/"
    else:
        ll = psp_pbesol
        base_folder = "../../../libraries-pbesol/"
        sssp_path = "./mix-sssp-prec-pbesol-lib-v2/library/"

    for ff in ll:
        f = ff.split(".")[0:-1]
        # element = f[0]
        ptyp = f[1]

        pcode = f[4]
        plib = f[5]
        pcomment = ".".join(f[6:])

        if ptyp == "nc" and plib == "dojo" and pcomment == "v0.4.1-std":
            lib_path = "nc-dojo-v0.4.1-std"
        elif ptyp == "nc" and plib == "dojo" and pcomment == "v0.4.1-str":
            lib_path = "nc-dojo-v0.4.1-str"
        elif ptyp == "nc" and plib == "dojo" and pcomment == "v0.5.0-std":
            lib_path = "nc-dojo-v0.5.0-std"
        elif ptyp == "nc" and plib == "spms":
            lib_path = "nc-spms-oncvpsp4"
        elif ptyp == "nc" and plib == "sg15":
            lib_path = "nc-sg15-oncvpsp4"
        elif ptyp == "paw" and plib == "psl" and "v1.0.0-high" in pcomment:
            lib_path = "paw-psl-v1.0.0-high"
        elif ptyp == "paw" and plib == "psl" and "v1.0.0-low" in pcomment:
            lib_path = "paw-psl-v1.0.0-low"
        elif ptyp == "paw" and plib == "psl" and "v0" in pcomment:
            lib_path = "paw-psl-v0.x"
        elif ptyp == "us" and plib == "psl" and "v1.0.0-high" in pcomment:
            lib_path = "us-psl-v1.0.0-high"
        elif ptyp == "us" and plib == "psl" and "v1.0.0-low" in pcomment:
            lib_path = "us-psl-v1.0.0-low"
        elif ptyp == "us" and plib == "psl" and "v0" in pcomment:
            lib_path = "us-psl-v0.x"
        elif ptyp == "us" and plib == "gbrv":
            lib_path = "us-gbrv-v1.x-upf2"
        elif ptyp == "paw" and plib == "jth" and pcomment == "v1.1-std":
            lib_path = "paw-jth-v1.1-std"
        elif ptyp == "paw" and plib == "jth" and pcomment == "v1.1-str":
            lib_path = "paw-jth-v1.1-str"
        elif ptyp == "paw" and plib == "wentzcovitch":
            lib_path = "paw-lanthanides-wentzcovitch"
        elif plib == "uni-marburg":
            lib_path = "paw-actinides-marburg"
        else:
            print(f"{ptyp}, {pcode}, {plib}, {pcomment}")
            raise ValueError(f"Unknown filename: '{ff}'")

        src = Path(f"{base_folder}") / f"{lib_path}" / ff
        if not src.exists():
            print(f"{ptyp}, {pcode}, {plib}, {pcomment}")
            raise FileNotFoundError(f"{src} not exist")

        dst = Path(sssp_path)
        dst.mkdir(exist_ok=True)

        dout = shutil.copy2(src, dst / f"{ff}")
        print(f"Success: copy to '{dout}'")


if __name__ == "__main__":
    extract("pbe")
    extract("pbesol")
