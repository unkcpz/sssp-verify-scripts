from pathlib import Path

from sssp_verify_scripts.controllers import ConvergenceEOSGroupSubmissionController

base_path = Path(__file__).resolve().parent.parent.parent / "libraries-pbe"
gbrv_lib_path = base_path / "us-gbrv-1.x-upf2"
jth_lib_path = base_path / "paw-jth-1.1-std"
dojo_lib_path = base_path / "nc-dojo-v0.4.1-std"
psl_paw_lib_path = base_path / "paw-psl-1.0.0-high"


# Hg, Ga, N, Cs, Mn
pseudos = [
    # GBRV-1.x
    #f"{gbrv_lib_path}/Hg.us.pbe.z_12.uspp.gbrv.v1.upf",
    #f"{gbrv_lib_path}/Ga.us.pbe.z_19.uspp.gbrv.v1.4.upf",
    #f"{gbrv_lib_path}/N.us.pbe.z_5.uspp.gbrv.v1.2.upf",
    f"{gbrv_lib_path}/Cs.us.pbe.z_9.uspp.gbrv.v1.upf",
    f"{gbrv_lib_path}/Mn.us.pbe.z_15.uspp.gbrv.v1.5.upf",
    ## JTH
    #f"{jth_lib_path}/Hg.paw.pbe.z_12.atompaw.jth.v1.1-std.upf",
    #f"{jth_lib_path}/Ga.paw.pbe.z_13.atompaw.jth.v1.1-std.upf",
    #f"{jth_lib_path}/N.paw.pbe.z_5.atompaw.jth.v1.1-std.upf",
    #f"{jth_lib_path}/Cs.paw.pbe.z_9.atompaw.jth.v1.1-std.upf",
    #f"{jth_lib_path}/Mn.paw.pbe.z_15.atompaw.jth.v1.1-std.upf",
    ## DOJO
    #f"{dojo_lib_path}/Hg.nc.pbe.z_20.oncvpsp3.dojo.v0.4.1-std.upf",
    #f"{dojo_lib_path}/Ga.nc.pbe.z_13.oncvpsp3.dojo.v0.4.1-std.upf",
    f"{dojo_lib_path}/N.nc.pbe.z_5.oncvpsp3.dojo.v0.4.1-std.upf",
    #f"{dojo_lib_path}/Cs.nc.pbe.z_9.oncvpsp3.dojo.v0.4.1-std.upf",
    #f"{dojo_lib_path}/Mn.nc.pbe.z_15.oncvpsp3.dojo.v0.4.1-std.upf",
    ## PSL-PAW
    #f"{psl_paw_lib_path}/Hg.paw.pbe.z_20.ld1.psl.v1.0.0-high.upf",
    #f"{psl_paw_lib_path}/Ga.paw.pbe.z_13.ld1.psl.v1.0.0-high.upf",
    #f"{psl_paw_lib_path}/N.paw.pbe.z_5.ld1.psl.v1.0.0-high.upf",
    #f"{psl_paw_lib_path}/Cs.paw.pbe.z_9.ld1.psl.v1.0.0-high.upf",
    #f"{psl_paw_lib_path}/Mn.paw.pbe.z_15.ld1.psl.v1.0.0-high.upf",
]

conf = 'DC'
protocol = 'balanced'

target_upf_lib = "experimental/upf/candidate"
comment = "EOS-convergence-protocol-compare"
computer = 'eiger-hq'
unit_num_cpus = 32
unit_memory_mb = 120000 # mb
unit_npool = 4

wavefunction_cutoff_list = [i for i in range(20, 201, 20)]

print("\n[bold purple]--- Launching Convergence EOS controller ---[/]\n")
ConvergenceEOSGroupSubmissionController(
    group_label="experimental/convergence/eos/balanced/dc",
    parent_group_label=target_upf_lib,
    pw_code=f"pw-7.0@{computer}",
    protocol=protocol,
    wavefunction_cutoff_list=wavefunction_cutoff_list,
    unit_num_cpus=unit_num_cpus,
    unit_memory_mb=unit_memory_mb,
    unit_npool=unit_npool,
    clean_workdir=True,
).submit_new_batch(verbose=True)
