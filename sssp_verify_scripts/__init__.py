import os

from aiida import orm

from aiida.plugins import WorkflowFactory, DataFactory
from aiida.engine import submit

UpfData = DataFactory('pseudo.upf')
VerificationWorkChain = WorkflowFactory('sssp_workflow.verification')

def submit_verification(pw_code, ph_code, upf, label, test_mode=False):
    if test_mode:
        protocol = 'test'
    else:
        protocol = 'theos'
    inputs = {
        'pw_code': pw_code,
        'ph_code': ph_code,
        'pseudo': upf,
        'label': orm.Str(label),
        'protocol': orm.Str(protocol),
        'properties_list': orm.List(list=[
            'delta_factor',
            # 'convergence:cohesive_energy',
            # 'convergence:phonon_frequencies',
            # 'convergence:pressure',
        ]),
        'options': orm.Dict(
                dict={
                    'resources': {
                        'num_machines': 1,
                        # 'num_cores': 8*2,
                        # 'memory_Mb': 1024*20*2,
                    },
                    'max_wallclock_seconds': 1200 * 4,
                    'withmpi': True,
                }),
        'parallelization': orm.Dict(dict={'npool': 4}),
        'clean_workdir_level': orm.Int(1),
    }

    node = submit(VerificationWorkChain, **inputs)
    return node

def verify_pseudos_in_folder(sssp_dir, element, pseudos_dict, pw_code, ph_code, test_mode=False):
    for key, value in pseudos_dict.items():
        pp_path = os.path.join(sssp_dir, element, key)
        label = value['label']
        if test_mode:
            label = f'{label}::T'
        with open(pp_path, 'rb') as stream:
            pseudo = UpfData(stream)

        node = submit_verification(pw_code, ph_code, pseudo, label, test_mode)
        node.description = label
        print(node, node.description)