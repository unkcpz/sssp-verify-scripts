import os
from pathlib import Path

from aiida import orm

Path.mkdir(Path.cwd().joinpath("oxygen_eos"), exist_ok=True)

pks = [
"588123", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.nc.z_6.oncvpsp3.dojo.v0.4.1-str (['measure.precision'])                       0
"588184", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.nc.z_6.oncvpsp4.sg15.v0 (['measure.precision'])                               0
"588258", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.paw.z_6.atompaw.jth.v1.1-std (['measure.precision'])                          0
"588344", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.paw.z_6.atompaw.jth.v1.1-str (['measure.precision'])                          0
"588403", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.paw.z_6.ld1.psl.v0.1 (['measure.precision'])                                  0
"588463", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.paw.z_6.ld1.psl.v1.0.0-high (['measure.precision'])                           0
"588531", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.paw.z_6.ld1.psl.v1.0.0-low (['measure.precision'])                            0
"588604", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.us.z_6.ld1.psl.v0.1 (['measure.precision'])                                   0
"588699", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.us.z_6.ld1.psl.v1.0.0-high (['measure.precision'])                            0
"588843", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.us.z_6.ld1.psl.v1.0.0-low (['measure.precision'])                             0
"588920", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.us.z_6.uspp.gbrv.v1.2 (['measure.precision'])                                 0
"589069", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.nc.z_6.oncvpsp3.dojo.v0.4.1-std (['measure.precision'])                       0
"589141", #  (Oxygen) (acwf at daint-mc-mrcloud-mem - default) O.nc.z_6.oncvpsp4.spms.v1 (['measure.precision'])                               0
]

for pk in pks:
    n = orm.load_node(pk)
    label = n.base.extras.all['label']
    label = label.split(' ')[-1]

    ecutwfc = int(n.inputs.measure.wavefunction_cutoff.value)
    ecutrho = int(n.inputs.measure.charge_density_cutoff.value)

    base_filename = f"{label}.ecut.wfc{ecutwfc}.rho{ecutrho}"

    command = f"cd oxygen_eos && aiida-sssp-workflow inspect {pk} -o {base_filename}"
    os.system(command) 