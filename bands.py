import aiida
from aiida.orm import load_node, Float, Bool
from aiida_sssp_workflow.workflows.convergence.bands import helper_bands_distence_difference

aiida.load_profile()

# nosym = False
band_parameters_a  = load_node(28746)
band_parameters_b  = load_node(28663)
band_structure_a = load_node(28744)
band_structure_b = load_node(28661)

res = helper_bands_distence_difference(
    band_structure_a, band_parameters_a,
    band_structure_b, band_parameters_b,
    Float(0.06122564129655),
    Float(5.0),
    Bool(True),
)

print(res.get_dict())

# nosym = True
band_parameters_a  = load_node(29426)
band_parameters_b  = load_node(29343)
band_structure_a = load_node(29424)
band_structure_b = load_node(29341)

res = helper_bands_distence_difference(
    band_structure_a, band_parameters_a,
    band_structure_b, band_parameters_b,
    Float(0.06122564129655),
    Float(5.0),
    Bool(True),
)

print(res.get_dict())