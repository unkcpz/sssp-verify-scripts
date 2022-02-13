import os

from aiida import orm

from aiida.plugins import WorkflowFactory, DataFactory
from aiida.engine import submit

UpfData = DataFactory('pseudo.upf')
VerificationWorkChain = WorkflowFactory('sssp_workflow.verification')

def submit_verification(pw_code, ph_code, upf, label):
    inputs = {
        'pw_code': pw_code,
        'ph_code': ph_code,
        'pseudo': upf,
        'label': orm.Str(label),
        'protocol_calculation': orm.Str('theos'),
        'protocol_criteria': orm.Str('precision'),
        'properties_list': orm.List(list=[
            'delta_factor',
            'convergence:cohesive_energy',
            'convergence:phonon_frequencies',
            'convergence:pressure',
        ]),
        'options': orm.Dict(
                dict={
                    'resources': {
                        'num_machines': 1,
                        # 'num_cores': 8*2,
                        # 'memory_Mb': 1024*20*2,
                    },
                    'max_wallclock_seconds': 1200,
                    'withmpi': True,
                }),
        'parallelization': orm.Dict(dict={'npool': 16}),
        'clean_workdir_level': orm.Int(1),
    }

    node = submit(VerificationWorkChain, **inputs)
    return node

def verify_pseudos_in_folder(sssp_dir, element, pseudos_dict, pw_code, ph_code):
    for key, value in pseudos_dict.items():
        pp_path = os.path.join(sssp_dir, element, key)
        label = value['label']

        with open(pp_path, 'rb') as stream:
            pseudo = UpfData(stream)

        node = submit_verification(pw_code, ph_code, pseudo, label)
        node.description = label
        print(node, node.description)