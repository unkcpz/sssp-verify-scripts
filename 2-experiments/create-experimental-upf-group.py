from pathlib import Path

from aiida import orm
from aiida.plugins import DataFactory

UpfData = DataFactory("pseudo.upf")

GROUP_NAME = "experimental/upf/candidate"

base_path = Path(__file__).resolve().parent / "libraries-pbe"
gbrv_lib_path = base_path / "us-gbrv-1.x-upf2"
jth_lib_path = base_path / "paw-jth-1.1-std"
dojo_lib_path = base_path / "nc-dojo-v0.4.1-std"
psl_paw_lib_path = base_path / "paw-psl-1.0.0-high"


# Hg, Ga, N, Cs, Mn
pseudos = [
    # GBRV-1.x
    f"{gbrv_lib_path}/Hg.us.pbe.z_12.uspp.gbrv.v1.upf",
    f"{gbrv_lib_path}/Ga.us.pbe.z_19.uspp.gbrv.v1.4.upf",
    f"{gbrv_lib_path}/N.us.pbe.z_5.uspp.gbrv.v1.2.upf",
    f"{gbrv_lib_path}/Cs.us.pbe.z_9.uspp.gbrv.v1.upf",
    f"{gbrv_lib_path}/Mn.us.pbe.z_15.uspp.gbrv.v1.5.upf",
    # JTH
    #f"{jth_lib_path}/Hg.paw.pbe.z_12.atompaw.jth.v1.1-std.upf",
    #f"{jth_lib_path}/Ga.paw.pbe.z_13.atompaw.jth.v1.1-std.upf",
    #f"{jth_lib_path}/N.paw.pbe.z_5.atompaw.jth.v1.1-std.upf",
    #f"{jth_lib_path}/Cs.paw.pbe.z_9.atompaw.jth.v1.1-std.upf",
    #f"{jth_lib_path}/Mn.paw.pbe.z_15.atompaw.jth.v1.1-std.upf",
    # DOJO
    f"{dojo_lib_path}/Hg.nc.pbe.z_20.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/Ga.nc.pbe.z_13.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/N.nc.pbe.z_5.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/Cs.nc.pbe.z_9.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/Mn.nc.pbe.z_15.oncvpsp3.dojo.v0.4.1-std.upf",
    # PSL-PAW
    f"{psl_paw_lib_path}/Hg.paw.pbe.z_20.ld1.psl.v1.0.0-high.upf",
    f"{psl_paw_lib_path}/Ga.paw.pbe.z_13.ld1.psl.v1.0.0-high.upf",
    f"{psl_paw_lib_path}/N.paw.pbe.z_5.ld1.psl.v1.0.0-high.upf",
    f"{psl_paw_lib_path}/Cs.paw.pbe.z_9.ld1.psl.v1.0.0-high.upf",
    f"{psl_paw_lib_path}/Mn.paw.pbe.z_15.ld1.psl.v1.0.0-high.upf",
]

group, _ = orm.Group.collection.get_or_create(label=GROUP_NAME)
group.description = """
This group contains the candidate UPF files for the experimental calculations.
Source data of the analysis in supplementary information of the paper.
It contains elements Hg, Ga, N, Cs, Mn from the following libraries:
- GBRV-1.x-upf2: Vanderbilt uspp code
- PAW-JTH-1.1-std: AtomPAW code
- NC-DOJO-v0.4.1-std: ONCVPSP code
- PAW-PSL-1.0.0-high: LD1 atomic code of QE
"""

for pseudo_path in pseudos:
    if not Path(pseudo_path).exists():
        raise FileNotFoundError(f"Pseudo file not found: {pseudo_path}")
    else:
        print(f"Pseudo file found: {pseudo_path}")

    pseudo = UpfData.get_or_create(pseudo_path).store()
    
    # check by query md5 in the group if pseudo is already in the group
    qb = orm.QueryBuilder()
    qb.append(
        orm.Group, filters={"label": GROUP_NAME}, tag="group"
    ).append(
        UpfData, filters={"extras.md5": pseudo.md5}, with_group="group"
    )
    
    if qb.count() != 0:
        print(f"pseudo {pseudo} is already in the group {GROUP_NAME}")
        continue

    pseudo.base.extras.set("md5", pseudo.md5)
    group.add_nodes(pseudo)
    print(f"Added {pseudo} to {GROUP_NAME}")