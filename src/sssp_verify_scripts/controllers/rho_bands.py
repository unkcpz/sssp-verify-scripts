from aiida_submission_controller import FromGroupSubmissionController

from aiida import orm
from aiida.engine import ProcessBuilder
from aiida_pseudo.data.pseudo import UpfData

from aiida_sssp_workflow.utils import extract_pseudo_info
from aiida_sssp_workflow.utils import get_default_dual
from aiida_sssp_workflow.utils.structure import UNARIE_CONFIGURATIONS
from aiida_sssp_workflow.workflows.convergence.bands import ConvergenceBandsWorkChain
from aiida_sssp_workflow.utils.pseudo import DualType, get_dual_type


class ConvergenceBandsGroupSubmissionController(FromGroupSubmissionController):
    """The submission controller for convergence bands group."""

    unique_extra_keys: tuple = ("md5",)
    parent_group_label: str
    group_label: str

    pw_code: str
    protocol: str
    element_wavefunction_cutoff_mapping: dict[str, float]
    unit_num_cpus: int
    unit_memory_mb: int
    unit_npool: int
    configuration: str = "DC"
    clean_workdir: bool

    def get_inputs_and_processclass_from_extras(self, extras_values):
        """Return the builder for the submission."""
        parent_node = self.get_parent_node_from_extras(extras_values)

        # the parent_node should be a pseudo node
        if not isinstance(parent_node, UpfData):
            raise ValueError(
                f"The parent node should be a UpfData node, but got {parent_node}"
            )

        if self.configuration not in UNARIE_CONFIGURATIONS:
            raise ValueError(
                f"Got {self.configuration}, the configuration shuold be on of {UNARIE_CONFIGURATIONS}."
            )

        pseudo = parent_node

        pp_info = extract_pseudo_info(pseudo.get_content())
        pseudo = parent_node
        element = pp_info.element
        dual = get_default_dual(pseudo)

        num_cpus = self.unit_num_cpus * 1
        memory_mb = self.unit_memory_mb * 1
        npool = self.unit_npool * 1

        dual_list = None
        match get_dual_type(pp_info.type, element):
            case DualType.NC:
                dual_list = [2.0, 2.5, 3.0, 3.5, 4.0]
            case DualType.AUGLOW:
                dual_list = [6.0, 6.5, 7.0, 7.5, 8.0]
            case DualType.AUGHIGH:
                dual_list = [8.0, 9.0, 10.0, 12.0, 16.0, 18.0]

                # atom_npool *= 2
                # atom_num_cpus *= 2
                # atom_memory_mb *= 2

        if dual_list is None:
            raise ValueError("dual_list can not be None")

        ecutwfc = self.element_wavefunction_cutoff_mapping[element]
        cutoff_list = [(ecutwfc, ecutwfc * dual) for dual in dual_list]

        builder: ProcessBuilder = ConvergenceBandsWorkChain.get_builder(
            pseudo=parent_node,
            protocol=self.protocol,
            cutoff_list=cutoff_list,
            configuration=self.configuration,
            code=orm.load_code(self.pw_code),
            parallelization={"npool": npool},
            mpi_options={
                "resources": {"num_machines": 1, "tot_num_mpiprocs": 48},
            },
            clean_workdir=self.clean_workdir,
        )

        return builder
