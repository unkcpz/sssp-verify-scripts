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
        
        builder: ProcessBuilder = ConvergenceEOSWorkChain.get_builder(
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

class TransferabilityEOSGroupSubmissionController(FromGroupSubmissionController):
    """The submission controller for transferability EOS group."""
    unique_extra_keys: tuple = ('md5',)
    parent_group_label: str
    group_label: str
    pw_code: str
    curate_type: str
    protocol: str
    ecutwfc: int | None = None
    cutoff_mapping: dict | None = None
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

        match get_dual_type(pp_info.type, element):
            case DualType.NC:
                dual = 4
            case DualType.AUGLOW:
                dual = 8
            case DualType.AUGHIGH:
                dual = 18
            case _:
                raise ValueError("Unknow condition to set dual")

        # the parent_node should be a pseudo node
        if not isinstance(parent_node, UpfData):
            raise ValueError(f"The parent node should be a UpfData node, but got {parent_node}")

        if self.ecutwfc is None and self.cutoff_mapping is None:
            raise ValueError("Cannot be both None")

        if self.ecutwfc is not None and self.cutoff_mapping is not None:
            raise ValueError("Can only set one of it to None")

        if self.ecutwfc is not None:
            cutoffs = (self.ecutwfc, self.ecutwfc * dual)

        # cutoff_mapping -> dict
        # {'Ag.paw......upf': {'md5': xxxx, 'cutoffs': (20, 20)}}
        elif self.cutoff_mapping is not None:
            for pp_name, info in self.cutoff_mapping.items():
                if info.get('md5') == extras_values[0]:
                    cutoffs = info.get('cutoffs')
                    print(f"[bold blue]Info:[/] Using cutoff={cutoffs} for {pp_name}.")

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


        builder: ProcessBuilder = TransferabilityEOSWorkChain.get_builder(
            code=code,
            pseudo=parent_node,
            protocol=self.protocol,
            curate_type=self.curate_type,
            cutoffs=cutoffs,
            parallelization={"npool": npool},
            mpi_options=mpi_options,
            clean_workdir=self.clean_workdir,
        )

        return builder
