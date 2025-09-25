from pathlib import Path
from typing import Tuple
from aiida_submission_controller import FromGroupSubmissionController

from aiida import orm
from aiida.engine import ProcessBuilder
from aiida_pseudo.data.pseudo import UpfData

from aiida_sssp_workflow.utils import DualType, extract_pseudo_info, get_dual_type
from aiida_sssp_workflow.utils.element import HIGH_DUAL_ELEMENTS
from aiida_sssp_workflow.utils.structure import UNARIE_CONFIGURATIONS
from aiida_sssp_workflow.workflows.convergence.eos import ConvergenceEOSWorkChain
from aiida_sssp_workflow.workflows.transferability.eos import TransferabilityEOSWorkChain

class ConvergenceEOSGroupSubmissionController(FromGroupSubmissionController):
    """The submission controller for convergence EOS group."""
    unique_extra_keys: tuple = ('md5',)
    parent_group_label: str
    group_label: str

    pw_code: str
    protocol: str
    element_wavefunction_cutoff_mapping: dict[str, float]
    unit_num_cpus: int
    unit_memory_mb: int
    unit_npool: int
    element_configuration_mapping: dict[str, str]
    clean_workdir: bool

    def get_inputs_and_processclass_from_extras(self, extras_values):
        """Return the builder for the submission."""
        parent_node = self.get_parent_node_from_extras(extras_values)

        # the parent_node should be a pseudo node
        if not isinstance(parent_node, UpfData):
            raise ValueError(f"The parent node should be a UpfData node, but got {parent_node}")

        pseudo = parent_node
        
        pp_info = extract_pseudo_info(pseudo.get_content())
        element = pp_info.element

        npool = self.unit_npool * 1

        dual_list = None
        match get_dual_type(pp_info.type, element):
            case DualType.NC:
                dual_list = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
            case DualType.AUGLOW:
                dual_list = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
            case DualType.AUGHIGH:
                dual_list = [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 16.0, 18.0]

                # atom_npool *= 2
                # atom_num_cpus *= 2
                # atom_memory_mb *= 2

        if dual_list is None:
            raise ValueError("dual_list can not be None")

        ecutwfc = self.element_wavefunction_cutoff_mapping[element]
        cutoff_list = [(ecutwfc, ecutwfc * dual) for dual in dual_list]

        configuration = self.element_configuration_mapping[element]
        
        builder: ProcessBuilder = ConvergenceEOSWorkChain.get_builder(
            pseudo=parent_node,
            protocol=self.protocol,
            cutoff_list=cutoff_list,
            configuration=configuration,
            code=orm.load_code(self.pw_code),
            parallelization={"npool": npool},
            mpi_options={
                'resources': {
                    'num_machines': 1,
                    'tot_num_mpiprocs': 48
                },
            },
            clean_workdir=self.clean_workdir,
        )

        return builder

