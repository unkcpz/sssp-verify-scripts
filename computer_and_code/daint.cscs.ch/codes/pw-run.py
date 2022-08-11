from aiida import orm, plugins
from aiida.engine import submit

code = orm.load_code('sarus-pw-6.8@daint-mc-mr0')
builder = code.get_builder()

structure = orm.load_node(120)
builder.structure = structure
pseudo_family = orm.load_group('SSSP/1.1/PBE/efficiency')
pseudos = pseudo_family.get_pseudos(structure=structure)
builder.pseudos = pseudos
parameters = {
  'CONTROL': {
    'calculation': 'scf',  # self-consistent field
  },
  'SYSTEM': {
    'ecutwfc': 30.,  # wave function cutoff in Ry
    'ecutrho': 240.,  # density cutoff in Ry
  },
}
builder.parameters = orm.Dict(dict=parameters)

KpointsData = plugins.DataFactory('array.kpoints')
kpoints = KpointsData()
kpoints.set_kpoints_mesh([4,4,4])
builder.kpoints = kpoints

builder.metadata.options.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 1}
builder.metadata.options.mpirun_extra_params = ['--mpi=pmi2']

calcjob_node = submit(builder)