import os

from aiida import orm

from aiida.plugins import WorkflowFactory, DataFactory
from aiida.engine import submit

from aiida_sssp_workflow.workflows.verifications import DEFAULT_PROPERTIES_LIST

UpfData = DataFactory('pseudo.upf')
VerificationWorkChain = WorkflowFactory('sssp_workflow.verification')

def run_verification(
    pw_code, ph_code, upf, properties_list=[], label=None, mode='test',
):
    """
    pw_code: code for pw.x calculation
    ph_code: code for ph.x calculation
    upf: upf file to verify
    properties_list: propertios to verified
    label: if None, label will parsed from filename
    mode: 
        test to run on localhost with test protocol
        precheck: precheck control protocol on convergence verification
        standard: running a production on eiger
    """
    clean_level = 1
    
    inputs = {
        "accuracy": {
            "protocol": orm.Str("test"),
            "cutoff_control": orm.Str("test"),
        },
        "convergence": {
            "protocol": orm.Str("test"),
            "cutoff_control": orm.Str("test"),
            "criteria": orm.Str("efficiency"),
            # "preset_ecutwfc": orm.Int(31),
        },
        "pw_code": pw_code,
        "ph_code": ph_code,
        "pseudo": upf,
        "properties_list": orm.List(list=properties_list),
        "label": orm.Str(label),
        "options": orm.Dict(
            dict={
                "resources": {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": 1,
                },
                "max_wallclock_seconds": 1800 * 3,
                "withmpi": True,
            }
        ),
        "parallelization": orm.Dict(dict={}),
        "clean_workdir_level": orm.Int(clean_level),
    }

    res, node = run_get_node(VerificationWorkChain, **inputs)
    return res, node

def verify_pseudos_in_folder(sssp_dir, element, pseudos_dict, pw_code, ph_code, test_mode=False):
    for key, value in pseudos_dict.items():
        pp_path = os.path.join(sssp_dir, element, key)
        dual = value['dual']
        label = value['label']
        if test_mode:
            label = f'{label}::T'
        with open(pp_path, 'rb') as stream:
            pseudo = UpfData(stream)

        node = submit_verification(pw_code, ph_code, pseudo, label, dual, test_mode)
        node.description = label
        print(node, node.description)
        
        
def verify_pseudos(sssp_path, element, pseudos_dict, pw_code, ph_code, test_mode=False):
    for key, value in pseudos_dict.items():
        pp_path = os.path.join(sssp_dir, element, key)
        dual = value['dual']
        label = value['label']
        if test_mode:
            label = f'{label}::T'
        with open(pp_path, 'rb') as stream:
            pseudo = UpfData(stream)

        node = submit_verification(pw_code, ph_code, pseudo, label, dual, test_mode)
        node.description = label
        print(node, node.description)