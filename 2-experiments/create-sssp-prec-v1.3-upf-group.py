from pathlib import Path

from aiida import orm
from aiida_pseudo.data.pseudo import UpfData

GROUP_NAME = "validate/upf/candidate/sssp-prec-v1.3"

base_path = Path(__file__).resolve().parent.parent / "libraries-pbe"
sssp_lib_path = base_path / "mix-sssp-prec-v1.3.0"

pseudos = [p.resolve() for p in sssp_lib_path.glob('*')]
print(len(pseudos))

group, _ = orm.Group.collection.get_or_create(label=GROUP_NAME)
group.description = """
This group contains the candidate UPF files from SSSP precision v1.3.0 for validate calculations.
Source data of the analysis in supplementary information of the paper.
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
