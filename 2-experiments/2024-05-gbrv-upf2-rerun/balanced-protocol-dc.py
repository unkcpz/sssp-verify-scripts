from sssp_verify_scripts.controllers import ConvergenceEOSGroupSubmissionController
from pydantic import ValidationError


conf = 'DC'
protocol = 'balanced'
lib_name = "us-gbrv-v1.x-upf2"

target_upf_lib = f"validate/upf/candidate/{lib_name}"
target_group_label=f"validate/{lib_name}/convergence/eos/{protocol}/{conf.lower()}"

comment = "EOS-convergence-protocol-compare"

computer = 'eiger-hq'
unit_num_cpus = 32
unit_memory_mb = 120000 # mb
unit_npool = 4

wavefunction_cutoff_list = [i for i in range(20, 201, 5)]

print("\nLaunching Convergence EOS controller ---\n")
try:
    controller = ConvergenceEOSGroupSubmissionController(
        group_label=target_group_label,
        parent_group_label=target_upf_lib,
        max_concurrent=5,
        pw_code=f"pw-7.0@{computer}",
        protocol=protocol,
        configuration=conf,
        wavefunction_cutoff_list=wavefunction_cutoff_list,
        unit_num_cpus=unit_num_cpus,
        unit_memory_mb=unit_memory_mb,
        unit_npool=unit_npool,
        clean_workdir=True,
    )
except ValidationError as exc:
    print(repr(exc.errors()))
else:
    controller.submit_new_batch(verbose=True)
