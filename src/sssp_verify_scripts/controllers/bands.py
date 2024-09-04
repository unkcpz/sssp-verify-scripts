from aiida_submission_controller import FromGroupSubmissionController

from aiida import orm
from aiida.engine import ProcessBuilder
from aiida_pseudo.data.pseudo import UpfData

from aiida_sssp_workflow.utils import get_default_dual
from aiida_sssp_workflow.utils.structure import UNARIE_CONFIGURATIONS
from aiida_sssp_workflow.workflows.convergence.bands import ConvergenceBandsWorkChain

class ConvergenceBandsGroupSubmissionController(FromGroupSubmissionController):
    """The submission controller for convergence bands group."""
    unique_extra_keys: tuple = ('md5',)
    parent_group_label: str
    group_label: str

    pw_code: str
    protocol: str
    wavefunction_cutoff_list: list
    unit_num_cpus: int
    unit_memory_mb: int
    unit_npool: int
    configuration: str = 'DC'
    clean_workdir: bool

    def get_inputs_and_processclass_from_extras(self, extras_values):
        """Return the builder for the submission."""
        parent_node = self.get_parent_node_from_extras(extras_values)

        # the parent_node should be a pseudo node
        if not isinstance(parent_node, UpfData):
            raise ValueError(f"The parent node should be a UpfData node, but got {parent_node}")

        if self.configuration not in UNARIE_CONFIGURATIONS:
            raise ValueError(f"Got {self.configuration}, the configuration shuold be on of {UNARIE_CONFIGURATIONS}.")

        pseudo = parent_node
        dual = get_default_dual(pseudo)

        num_cpus = self.unit_num_cpus * 1
        memory_mb = self.unit_memory_mb * 1
        npool = self.unit_npool * 1

        dual = 12
        # cutoff_list = [(ecutwfc, ecutwfc * dual) for ecutwfc in self.wavefunction_cutoff_list]
        cutoff_list = [(ecutwfc, ecutwfc * dual) for ecutwfc in self.wavefunction_cutoff_list[:-1]] + [(200, 3600)]
        
        builder: ProcessBuilder = ConvergenceBandsWorkChain.get_builder(
            pseudo=parent_node,
            protocol=self.protocol,
            cutoff_list=cutoff_list,
            configuration=self.configuration,
            code=orm.load_code(self.pw_code),
            parallelization={"npool": npool},
            mpi_options={
                'resources': {
                    'num_cpus': num_cpus,
                    'memory_mb': memory_mb,
                },
            },
            clean_workdir=self.clean_workdir,
        )

        return builder
