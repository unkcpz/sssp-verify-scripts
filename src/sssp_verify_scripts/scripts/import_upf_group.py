from pathlib import Path
import argparse

import aiida
from aiida import orm
from aiida_pseudo.data.pseudo import UpfData

def get_lib_info(lib_name: str, base_path: Path) -> tuple:
    if lib_name == 'nc-dojo-v0.4.1-std':
        description = """This group contains the candidate UPF files from DOJO std 0.4.1 library for validate calculations.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'nc-dojo-v0.4.1-str':
        description = """This group contains the candidate UPF files from DOJO str 0.4.1 library for validate calculations.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'nc-dojo-v0.5.0-std':
        description = """This group contains the candidate UPF files from DOJO std 0.5.0 library for validate calculations.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'nc-psl-v1.0.0':
        description = """This group contains the candidate UPF files from PSLibrary NC v1.0.0 library for validate calculations.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'nc-sg15-oncvpsp4':
        description = """This group contains the candidate UPF files by ONCVPSP from SG15 library for validate calculations.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'nc-spms-oncvpsp4':
        description = """This group contains the candidate UPF files by ONCVPSP from SPMS library for validate calculations.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'paw-jth-v1.1-std':
        description = """This group contains the candidate UPF files from JTH std v1.1 library for validate calculations.
PPs are regenerated with PP pseudo orbital label for every channel.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'paw-jth-v1.1-str':
        description = """This group contains the candidate UPF files from JTH str v1.1 library for validate calculations.
PPs are regenerated with PP pseudo orbital label for every channel.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'paw-psl-v1.0.0-high':
        description = """This group contains the candidate UPF files from PSLibrary PAW v1.0.0 high library for validate calculations.
PPs are regenerated with inputs.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'paw-psl-v1.0.0-low':
        description = """This group contains the candidate UPF files from PSLibrary PAW v1.0.0 low library for validate calculations.
PPs are regenerated with inputs.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'paw-lanthanides-wentzcovitch':
        description = """This group contains the candidate UPF files by ATOMPAW from Wentzcovitch library low library for validate calculations.
PPs are regenerated with inputs.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'paw-actinides-marburg':
        description = """This group contains the candidate UPF files by ATOMPAW from Marburg Uni library low library for validate calculations.
PPs are regenerated with inputs.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'paw-psl-v0.x':
        description = """This group contains the candidate UPF files from PSLibrary PAW v0.x library for validate calculations.
PPs are regenerated with inputs.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'us-psl-v1.0.0-high':
        description = """This group contains the candidate UPF files from PSLibrary US v1.0.0 high library for validate calculations.
PPs are regenerated with inputs.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'us-psl-v1.0.0-low':
        description = """This group contains the candidate UPF files from PSLibrary US v1.0.0 low library for validate calculations.
PPs are regenerated with inputs.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'us-psl-v0.x':
        description = """This group contains the candidate UPF files from PSLibrary US v0.x library for validate calculations.
PPs are regenerated with inputs.
"""
        group_name = f"validate/upf/candidate/{lib_name}"
        folder_path = base_path / "libraries-pbe" / f"{lib_name}"
    elif lib_name == 'us-gbrv-v1.x-upf2':
        description = """This group contains the candidate UPF files from GBRV library for validate calculations.
The PPs were regenerated with to UPF2 format with PP_QIJL,
Source data of the analysis in supplementary information of the paper.
"""
        group_name = "validate/upf/candidate/us-gbrv-v1.x-upf2"
        folder_path = base_path / "libraries-pbe" / "us-gbrv-v1.x-upf2"
    elif lib_name == 'mix-sssp-prec-v1.3.0':
        # FIXME: and rerun, use latest USPP PP and JTH PP.
        description = """This group contains the candidate UPF files from SSSP precision v1.3.0 for validate calculations.
Source data of the analysis in supplementary information of the paper.
"""
        group_name = "validate/upf/candidate/sssp-prec-v1.3"
        folder_path = base_path / "libraries-pbe" / "mix-sssp-prec-v1.3.0"
    elif lib_name == 'high-dual-elements':
        description = """High dual elements for test purpose, Fe, O, Hf
"""
        group_name = "validate/upf/candidate/high-dual-elements"
        folder_path = base_path / "libraries-pbe" / "high-dual-elements"
    else:
        raise ValueError(f'{lib_name} not found')

    return group_name, folder_path, description

def import_upf(group_name: str, folder: Path, description: str, dry: bool):
    pseudos = [p.resolve() for p in folder.glob('*.upf')]
    print(f'There are {len(pseudos)} UPF files in this folder')
    
    if dry:
        print(f'Importing to {group_name} from source folder {folder.resolve()}.')
        print(description)
        print('The pseudos that will be imported are:')
        print('\n'.join([p.name for p in pseudos]))

        return 

    group, _ = orm.Group.collection.get_or_create(label=group_name)
    group.description = description
    
    for pseudo_path in pseudos:
        if not Path(pseudo_path).exists():
            raise FileNotFoundError(f"Pseudo file not found: {pseudo_path}")
        else:
            print(f"Pseudo file found: {pseudo_path}")
    
        pseudo = UpfData.get_or_create(pseudo_path).store()
        
        # check by query md5 in the group if pseudo is already in the group
        qb = orm.QueryBuilder()
        qb.append(
            orm.Group, filters={"label": group_name}, tag="group"
        ).append(
            UpfData, filters={"extras.md5": pseudo.md5}, with_group="group"
        )
        
        if qb.count() != 0:
            print(f"pseudo {pseudo} is already in the group {group_name}")
            continue
    
        pseudo.base.extras.set("md5", pseudo.md5)
        group.add_nodes(pseudo)
        print(f"Added {pseudo} to {group_name}")

def main():
    aiida.load_profile()

    parser = argparse.ArgumentParser()
    parser.add_argument('--lib-name', help="The lib name.")
    parser.add_argument('--dry', action='store_true', help="dry run")
    parser.add_argument('--base-folder', default=".", help="The base path where all libraries exist.")

    args = parser.parse_args()

    base_folder = Path(args.base_folder).resolve()
    print(f"The base path is {base_folder}")

    group_name, folder, description = get_lib_info(args.lib_name, base_folder)
    import_upf(group_name=group_name, folder=folder, description=description, dry=args.dry)

if __name__ == "__main__":
    main()
