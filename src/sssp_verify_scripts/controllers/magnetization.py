from aiida_submission_controller import FromGroupSubmissionController

from aiida import orm
from aiida.engine import ProcessBuilder
from aiida_pseudo.data.pseudo import UpfData

from aiida_sssp_workflow.utils import DualType, extract_pseudo_info, get_dual_type
from aiida_sssp_workflow.workflows.list_run.magnetization import MagnetizationChangeWorkChain

class MagnetizationGroupSubmissionController(FromGroupSubmissionController):
    """The submission controller for transferability EOS group."""
    unique_extra_keys: tuple = ('md5',)
    parent_group_label: str
    group_label: str
    configuration: str
    pw_code: str
    protocol: str
    scale_list: list
    ecutwfc: int
    unit_num_cpus: int
    unit_memory_mb: int
    unit_npool: int
    clean_workdir: bool

    def get_inputs_and_processclass_from_extras(self, extras_values):
        """Return the builder for the submission."""
        parent_node = self.get_parent_node_from_extras(extras_values)

        pseudo: UpfData = parent_node
        
        pp_info = extract_pseudo_info(pseudo.get_content())
        element = pp_info.element

        # the parent_node should be a pseudo node
        if not isinstance(parent_node, UpfData):
            raise ValueError(f"The parent node should be a UpfData node, but got {parent_node}")

        num_cpus = self.unit_num_cpus * 1
        memory_mb = self.unit_memory_mb * 1
        npool = self.unit_npool * 1

        code = orm.load_code(self.pw_code)
        if 'hq' in code.computer.label:
            mpi_options = {
                'resources': {
                    'num_cpus': num_cpus,
                    'memory_mb': memory_mb,
                },
            }
        else:
            mpi_options={
                'resources': {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": num_cpus,
                },
            }

        match get_dual_type(pp_info.type, element):
            case DualType.NC:
                dual = 4
            case DualType.AUGLOW:
                dual = 8
            case DualType.AUGHIGH:
                dual = 16

        cutoffs = (self.ecutwfc, self.ecutwfc * dual)

        builder: ProcessBuilder = MagnetizationChangeWorkChain.get_builder(
            code=code,
            pseudo=parent_node,
            protocol=self.protocol,
            scale_list=self.scale_list,
            cutoffs=cutoffs,
            configuration=self.configuration,
            parallelization={"npool": npool},
            mpi_options=mpi_options,
            clean_workdir=self.clean_workdir,
        )

        return builder
