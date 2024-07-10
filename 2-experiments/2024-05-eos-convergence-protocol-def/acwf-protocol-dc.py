from sssp_verify_scripts.controllers import ConvergenceEOSGroupSubmissionController
from pydantic import ValidationError


conf = 'DC'
protocol = 'acwf'
mode = "validate"

if mode == "experimental":
    # Experiments on the candidate PPs: 4 * 5
    target_upf_lib = "experimental/upf/candidate"
    target_group_label=f"experimental/convergence/eos/{protocol}/{conf.lower()}"
elif mode == "validate":    
    # Validate on the candidate PPs: SSSP-prec-v1.3 recollected
    target_upf_lib = "validate/upf/candidate/sssp-prec-v1.3"
    target_group_label=f"validate/sssp-prec-v1.3/convergence/eos/{protocol}/{conf.lower()}"
else:
    raise ValueError(f"Unknown mode {mode}")

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
