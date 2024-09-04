from aiida_submission_controller import FromGroupSubmissionController

from aiida import orm
from aiida.engine import ProcessBuilder
from aiida_pseudo.data.pseudo import UpfData

from aiida_sssp_workflow.utils import extract_pseudo_info
from aiida_sssp_workflow.utils.element import HIGH_DUAL_ELEMENTS
from aiida_sssp_workflow.utils.structure import UNARIE_CONFIGURATIONS
from aiida_sssp_workflow.workflows.convergence.cohesive_energy import ConvergenceCohesiveEnergyWorkChain
from aiida_sssp_workflow.utils.pseudo import DualType, get_dual_type

class ConvergenceCohesiveEnergyGroupSubmissionController(FromGroupSubmissionController):
    """The submission controller for convergence EOS group."""
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
        
        pp_info = extract_pseudo_info(pseudo.get_content())
        element = pp_info.element

        num_cpus = self.unit_num_cpus * 1
        memory_mb = self.unit_memory_mb * 1
        npool = self.unit_npool * 1

        atom_num_cpus = num_cpus
        atom_memory_mb = memory_mb
        atom_npool = 1  # For isolate atom always use 1 npool

        match get_dual_type(pp_info.type, element):
            case DualType.NC:
                dual = 4
            case DualType.AUGLOW:
                dual = 8
            case DualType.AUGHIGH:
                dual = 18

                # atom_npool *= 2
                # atom_num_cpus *= 2
                # atom_memory_mb *= 2


        dual = 12
        # cutoff_list = [(ecutwfc, ecutwfc * dual) for ecutwfc in self.wavefunction_cutoff_list]
        cutoff_list = [(ecutwfc, ecutwfc * dual) for ecutwfc in self.wavefunction_cutoff_list[:-1]] + [(200, 3600)]
        
        builder: ProcessBuilder = ConvergenceCohesiveEnergyWorkChain.get_builder(
            pseudo=parent_node,
            protocol=self.protocol,
            cutoff_list=cutoff_list,
            configuration=self.configuration,
            code=orm.load_code(self.pw_code),
            bulk_parallelization={"npool": npool},
            bulk_mpi_options={
                'resources': {
                    'num_cpus': num_cpus,
                    'memory_mb': memory_mb,
                },
            },
            atom_parallelization={"npool": atom_npool},
            atom_mpi_options={
                'resources': {
                    'num_cpus': atom_num_cpus,
                    'memory_mb': atom_memory_mb,
                },
            },
            clean_workdir=self.clean_workdir,
        )

        return builder
