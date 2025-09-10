from aiida_submission_controller import FromGroupSubmissionController

from aiida import orm
from aiida.engine import ProcessBuilder
from aiida_pseudo.data.pseudo import UpfData

from aiida_sssp_workflow.utils import extract_pseudo_info
from aiida_sssp_workflow.utils.element import HIGH_DUAL_ELEMENTS
from aiida_sssp_workflow.utils.structure import UNARIE_CONFIGURATIONS
from aiida_sssp_workflow.workflows.convergence.phonon_frequencies import ConvergencePhononFrequenciesWorkChain

class ConvergencePhononFrequenciesGroupSubmissionController(FromGroupSubmissionController):
    """The submission controller for convergence EOS group."""
    unique_extra_keys: tuple = ('md5',)
    parent_group_label: str
    group_label: str

    pw_code: str
    ph_code: str
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
        
        pp_info = extract_pseudo_info(pseudo.get_content())
        element = pp_info.element

        if pp_info.type == 'nc':
            dual = 4
        else:
            dual = 8

        if element in HIGH_DUAL_ELEMENTS and pp_info.type != 'nc':
            dual = 18

            num_cpus = self.unit_num_cpus * 2
            memory_mb = self.unit_memory_mb * 2
            npool = self.unit_npool * 2
        else:
            num_cpus = self.unit_num_cpus * 1
            memory_mb = self.unit_memory_mb * 1
            npool = self.unit_npool * 1

        cutoff_list = [(ecutwfc, ecutwfc * dual) for ecutwfc in self.wavefunction_cutoff_list]
        
        builder: ProcessBuilder = ConvergencePhononFrequenciesWorkChain.get_builder(
            pseudo=parent_node,
            protocol=self.protocol,
            cutoff_list=cutoff_list,
            configuration=self.configuration,
            pw_code=orm.load_code(self.pw_code),
            ph_code=orm.load_code(self.ph_code),
            pw_parallelization={"npool": npool},
            pw_mpi_options={
                'resources': {
                    'num_cpus': num_cpus,
                    'memory_mb': memory_mb,
                },
            },
            ph_settings = {"CMDLINE": ["-npool", f"{npool}"]},
            ph_mpi_options={
                'resources': {
                    'num_cpus': num_cpus,
                    'memory_mb': memory_mb,
                },
            },
            clean_workdir=self.clean_workdir,
        )

        return builder
