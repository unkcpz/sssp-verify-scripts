import os
from pathlib import Path

from aiida import orm

Path.mkdir(Path.cwd().joinpath("nitrigen_eos"), exist_ok=True)

pks = [
"584942", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.paw.z_5.atompaw.jth.v1.1-std (['measure.precision'])                        0
"584997", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.nc.z_5.oncvpsp3.dojo.v0.4.1-std (['measure.precision'])                     0
"587104", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.nc.z_5.oncvpsp4.sg15.v0 (['measure.precision'])                             0
"587119", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.nc.z_5.oncvpsp4.spms.v1 (['measure.precision'])                             0
"587187", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.paw.z_5.atompaw.jth.v1.1-str (['measure.precision'])                        0
"587273", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.paw.z_5.ld1.psl.v0.1 (['measure.precision'])                                0
"587327", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.paw.z_5.ld1.psl.v1.0.0-high (['measure.precision'])                         0
"587363", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.us.z_5.ld1.psl.v0.1 (['measure.precision'])                                 0
"587476", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.us.z_5.ld1.psl.v1.0.0-high (['measure.precision'])                          0
"587595", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.us.z_5.ld1.theose.v0 (['measure.precision'])                                0
"587700", #  (Nitrigen) (acwf at daint-mc-mrcloud-mem - default) N.us.z_5.uspp.gbrv.v1.2 (['measure.precision'])                               0
]

for pk in pks:
    n = orm.load_node(pk)
    label = n.base.extras.all['label']
    label = label.split(' ')[-1]

    ecutwfc = int(n.inputs.measure.wavefunction_cutoff.value)
    ecutrho = int(n.inputs.measure.charge_density_cutoff.value)

    base_filename = f"{label}.ecut.wfc{ecutwfc}.rho{ecutrho}"

    command = f"cd nitrigen_eos && aiida-sssp-workflow inspect {pk} -o {base_filename}"
    os.system(command) 