import numpy as np

import aiida
from tqdm import tqdm

from aiida import orm
from aiida_sssp_workflow.workflows.convergence.report import ConvergenceReport
from aiida_sssp_workflow.calculations.calculate_metric import rel_errors_vec_length

class NuCutoffCalculator:
    
    criteria: float

    def __init__(self, criteria):
        self.criteria = criteria

    def extract(self, sample_uuid, ref_uuid) -> float:
        sample_node = orm.load_node(sample_uuid)
        ref_node = orm.load_node(ref_uuid)

        V0, B0, B1 = sample_node.outputs.output_parameters.get_dict()["birch_murnaghan_results"]
        ref_V0, ref_B0, ref_B1 = ref_node.outputs.output_parameters.get_dict()["birch_murnaghan_results"]

        nu = rel_errors_vec_length(ref_V0, ref_B0, ref_B1, V0, B0, B1)

        return nu

class EpsilonCutoffCalculator:
    """Epsilon metric from nat.rev.phys, the average energies of all(usually 7) EOS volume points"""

    criteria: float

    def __init__(self, criteria):
        self.criteria = criteria

    def extract(self, sample_uuid, ref_uuid) -> float:
        sample_node = orm.load_node(sample_uuid)
        ref_node = orm.load_node(ref_uuid)

        arr_sample = np.array(sample_node.outputs.eos.output_volume_energy.get_dict()["energies"])
        arr_ref = np.array(ref_node.outputs.eos.output_volume_energy.get_dict()["energies"])

        avg_sample = np.average(arr_sample)
        avg_ref = np.average(arr_ref)

        # eq.6 of nat.rev.phys
        A = np.sum(np.square(arr_sample - arr_ref))
        B = np.sum(np.square(arr_sample - avg_sample))
        C = np.sum(np.square(arr_ref - avg_ref))

        epsilon = np.sqrt(A / (np.sqrt(B * C)))

        return epsilon

    
def compute_recommended_cutoffs(node: orm.WorkChainNode, calculator) -> tuple:
    report = ConvergenceReport.construct(**node.outputs.convergence_report)

    ref = report.reference
    for p in reversed(report.convergence_list):
        if p.exit_status != 0:
            continue
        
        try:
            value = calculator.extract(p.uuid, ref.uuid)
        except Exception:
            continue

        if value > calculator.criteria:
            v = round(value, 4)
            return (p.wavefunction_cutoff, p.charge_density_cutoff, v)

    # Which means the first cutoff already converged.
    return (p.wavefunction_cutoff, p.charge_density_cutoff, -1.0)
        


def main():
    balanced_group: orm.Group = orm.Group.collection.get(label="validate/sssp-prec-v1.3/convergence/eos/balanced/dc")
    acwf_group: orm.Group = orm.Group.collection.get(label="validate/sssp-prec-v1.3/convergence/eos/acwf/dc")

    #calculator = NuCutoffCalculator(criteria=0.02)
    calculator = EpsilonCutoffCalculator(criteria=0.02)

    balanced_data = {}
    acwf_data = {}

    for gp, data in zip([balanced_group, acwf_group], [balanced_data, acwf_data]):
        for node in tqdm(gp.nodes):
            element = node.inputs.pseudo.element

            md5 = node.base.extras.get('md5', None)
            if md5 is None:
                raise ValueError(f'Node uuid={node.uuid} has no md5 extras.')

            if node.exit_status != 0:
                print(f"The workflow pk={node.pk} not finished okay, pseudo is {node.inputs.pseudo.filename}")
                data[md5] = (None, node.pk, element)
                continue

            wavefunction_cutoff, _, value = compute_recommended_cutoffs(node, calculator) 
            data[md5] = (wavefunction_cutoff, node.pk, element, value)

    for k in balanced_data.keys():
        if not (k in balanced_data and k in acwf_data):
            continue

        try:
            balanced_cutoff, balanced_pk, element, balanced_value = balanced_data.get(k)
            acwf_cutoff, acwf_pk, element, acwf_value = acwf_data.get(k)
        except:
            balanced_cutoff = None
            acwf_cutoff = None
            balanced_pk = None
            acwf_pk = None
            balanced_value = None
            acwf_value = None

        print(f"(balenced, acwf) -> cutoffs = ({balanced_cutoff}, {acwf_cutoff}) -> pks = ({balanced_pk}, {acwf_pk}), -> values = ({balanced_value}, {acwf_value}), element = {element}")

if __name__ == '__main__':
    main()
