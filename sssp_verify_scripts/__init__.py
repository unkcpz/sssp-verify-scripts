import os

from aiida import orm

from aiida.plugins import WorkflowFactory
from aiida.engine import submit, run_get_node

VerificationWorkChain = WorkflowFactory('sssp_workflow.verification')

def run_verification(
    pseudo,
    pw_code, 
    ph_code, 
    protocol, 
    cutoff_control, 
    criteria,
    options,
    parallization,
    properties_list,
    label,
    extra_desc,
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
    if cutoff_control.value == 'precheck':
        clean_level = 1  # not clean bands workchain and phonon frequencies
    elif cutoff_control.value == 'standard':
        clean_level = 2 # phonon frequencies workchan clean only finished okay.
    else:
        clean_level = 0
    
    inputs = {
        "accuracy": {
            "protocol": protocol,
            "cutoff_control": cutoff_control,
        },
        "convergence": {
            "protocol": protocol,
            "cutoff_control": cutoff_control,
            "criteria": criteria,
            # "preset_ecutwfc": orm.Int(31),
        },
        "pw_code": pw_code,
        "ph_code": ph_code,
        "pseudo": pseudo,
        "label": label,
        "properties_list": properties_list,
        "options": options,
        "parallelization": parallization,
        "clean_workdir_level": orm.Int(clean_level),  
    }
    
    if cutoff_control.value == 'test':
        _, node = run_get_node(VerificationWorkChain, **inputs)
    else:
        node = submit(VerificationWorkChain, **inputs)
        
    node.description = f'{label.value} ({extra_desc})'
        
    return node
