#!/usr/bin/env python
import json
import matplotlib.pyplot as plt
import os
import sys
import copy

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})


CAPPROPS = {
    'linewidth': 1,
    'color': 'tab:green'
}
BOXPROPS = {
    'linewidth': 1,
    'color': 'tab:green'
}
WHISKERPROPS = {
    'linewidth': 1,
    'color': 'tab:green'
}
MEDIANPROPS = {
    'linewidth': 1,
    'color': 'tab:blue'
}
FLIERPROPS = {
    'marker': '.',
    'linestyle': '',
    'markerfacecolor': 'tab:grey',
    'markeredgewidth': 0,
    'markersize': 3
}
# %%
def get_recommended_ecutwfc(element_data):
    """
    Get the recommended ecutwfc from the metadata
    """
    return element_data['cutoff_wfc']

def get_recommended_ecutrho(element_data):
    """
    Get the recommended ecutrho from the metadata
    """
    return element_data['cutoff_rho']

quantity_for_comparison_map = {
    "(Ry) wavefunction cutoff": get_recommended_ecutwfc,
    "(Ry) charge density cutoff": get_recommended_ecutrho,
}

xlims = {
    "(Ry) wavefunction cutoff": [0, 200.0],
    "(Ry) charge density cutoff": [0, 1600.0],
}

quantity_names = ["(Ry) wavefunction cutoff", "(Ry) charge density cutoff"]

QUANTITY_FANCY_NAMES = {
    'ecutwfc': "Wavefunction cutoff (Ry)",
    'ecutrho': "Charge density cutoff (Ry)",
}

DATA_FOLDER = "../lib-data/convergence"
with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
    labels_data = json.load(fhandle)
lib_labels = list(labels_data['methods-main'].keys())[::-1] # invert order because they are plot bottom to top, so we keep them alphabetical


def generate_box_plt(file_name, material_set_label, file_suffix, only_must_have_elements=None, keep_only_libs=None):
    """
    Generate the box plot
    """
    # Double check that there are no mistakes
    if keep_only_libs is not None:
        for lib in keep_only_libs:
            if lib not in lib_labels:
                raise ValueError(f"Asking to keep lib '{lib}' but it does not exist a lib with such a label")

    # Get the filtered list
    used_lib_labels = []
    for lib_label in lib_labels:
        if keep_only_libs is not None and lib_label not in keep_only_libs:
            print(f">> Skipping lib {lib_label}")
            continue
        used_lib_labels.append(lib_label)

    all_data = {}
    print()
    print('#####################################################################')
    print('#      Statistics on the number of elements for each lib            #')
    print('#####################################################################')
    for quantity_name in quantity_names:
        out_data = {}
        for lib_label in used_lib_labels:
            plugin_values = []
            plugin_big = 0
            plugin_small = 0
            out_data[lib_label] = {}

            with open(os.path.join(DATA_FOLDER, labels_data['methods-main'][lib_label]['metadata'])) as fhandle:
                pp_lib_data = json.load(fhandle)

            data = pp_lib_data

            plot_systems = set(
                el for el in data.keys()
                if data[el] is not None
            )
                           
            missing = []
            new_plot_systems = set()
            for element in only_must_have_elements:
                expected_key = f"{element}"
                if expected_key not in plot_systems:
                    missing.append(expected_key)
                else:
                    new_plot_systems.add(expected_key)

            plot_systems = new_plot_systems

            for element in plot_systems:

                element_data = data[element]

                try:
                    quantity_value = quantity_for_comparison_map[quantity_name](
                        element_data,
                    )
                except KeyError:
                    print(f"Skipping {element} because no data for {quantity_name}")
                    continue

                plugin_values.append(quantity_value)
                if quantity_value < xlims[quantity_name][0]:
                    plugin_big = plugin_big+1  
                if quantity_value > xlims[quantity_name][1]:
                    plugin_small = plugin_small+1

            out_data[lib_label]['values'] = plugin_values
            out_data[lib_label]['big'] = plugin_big
            out_data[lib_label]['small'] = plugin_small
        all_data[quantity_name] = out_data
    # %%
    # Set up the plot axes for each quantity
    print()
    print('#####################################################################')
    print('#   Description of the outliers (out of picture) for each lib       #')
    print('#####################################################################')
    for yy in quantity_names:
        for lib_label in used_lib_labels:
            if all_data[yy][lib_label]['small'] != 0:
                print(f'small {yy} {lib_label}', all_data[yy][lib_label]['small'])
    for yy in quantity_names:
        for lib_label in used_lib_labels:
            if all_data[yy][lib_label]['big'] != 0:
                print(f'big {yy} {lib_label}', all_data[yy][lib_label]['big'])
    
    n_quantities = len(quantity_names)
    # + 2 is to keep space for header and footer
    fig_height = max((len(used_lib_labels) + 2) * 0.45, 3) * 0.7 # Set min size of 3; rescaling factor to make it less wide
    fig, axes = plt.subplots(1, n_quantities, dpi=300, figsize=(4 * n_quantities, fig_height), sharey=True)
    axes = axes.flatten()

    for quantity_name, ax in zip(quantity_names, axes):
        quantity_values = [all_data[quantity_name][plugin_name]['values'] for plugin_name in used_lib_labels]
    
        fancy_labels = []
        for label in used_lib_labels:
            lib_label, sep, rest = label.partition('@')
            # Use bold code name, and replace pipe as with OT1 font enc it is replaced by a dash
            fancy_labels.append(rf'\textbf{{{lib_label}}}{sep}{rest}'.replace('|', '$|$'))

        ax.boxplot(
            quantity_values,
            flierprops=FLIERPROPS,
            boxprops=BOXPROPS,
            whiskerprops=WHISKERPROPS,
            medianprops=MEDIANPROPS,
            capprops=CAPPROPS,
            vert=False,
            labels=fancy_labels
        )
        ax.set_xlabel(rf"{quantity_name}",fontsize=14)
        ax.tick_params(axis='x', which='major', labelsize=11)
        ax.tick_params(axis='y', which='major', labelsize=9)
        ax.set_xlim(xlims[quantity_name][0],xlims[quantity_name][1])

    fig.suptitle(f"Materials set: {material_set_label}",fontsize=14, y=0.98)
    #fig.tight_layout()
    fig.subplots_adjust(left=0.25, right=0.99, top=1., bottom=(1.8/(len(used_lib_labels) + 2)), wspace=0.05)
    def make_space_above(axes, topmargin=1):
        """ increase figure size to make topmargin (in inches) space for 
            titles, without changing the axes sizes"""
        fig = axes.flatten()[0].figure
        s = fig.subplotpars
        w, h = fig.get_size_inches()

        figh = h - (1-s.top)*h  + topmargin
        fig.subplots_adjust(bottom=s.bottom*h/figh, top=1-topmargin/figh)
        fig.set_figheight(figh)
    make_space_above(axes, topmargin=0.35)
    fig.savefig(f'{file_name}{file_suffix}.pdf')


if __name__ == "__main__":
    import ase.data

    try:
        elements = sys.argv[1]
    except IndexError:
        print(
            "Pass as second parameter 'all' (atomic number 1-96), 'up-to-Bi-no-lanthanides' (1-56,71-83), "
            "'delta-set' (1-56,71-84+86), 'no-actinides' (1-86), 'only-actinides' (84-96), 'only-lanthanides'(57-71)."
        )
        sys.exit(1)

    only_must_have_elements = None
    if elements == 'all':
        chemical_numbers = list(range(1, 96+1))
        material_set_label = "Z=1-96"
    elif elements == 'delta-set':
        chemical_numbers = list(range(1, 56+1)) +  list(range(71, 84+1)) + [86]
        material_set_label = "Delta set (Science 2016)"
    elif elements == 'up-to-Bi-no-lanthanides':
        chemical_numbers = list(range(1, 56+1)) +  list(range(72, 83+1))
        material_set_label = "Z=1-56,72-83 (H to Ba, and Hf to Bi)"
    elif elements == 'no-actinides':
        chemical_numbers = list(range(1, 88+1))
        material_set_label = "Z=1-88"
    elif elements == 'only-actinides':
        chemical_numbers = list(range(84, 96+1))
        material_set_label = "Z=84-96 (Po to Cm)"
    elif elements == 'only-lanthanides':
        chemical_numbers = list(range(57, 71+1))
        material_set_label = "Z=57-71 (lanthanides: La to Lu)"
    else:
        print(
            "Pass as second parameter 'all' (atomic number 1-96), 'up-to-Bi-no-lanthanides' (1-56,71-83), "
            "'delta-set' (1-56,71-84+86), 'no-actinides' (1-86), 'only-actinides' (84-96), 'only-lanthanides'(57-71)."
        )
        sys.exit(1)
    # VASP  chemical_numbers.remove(1) #H
    #       chemical_numbers.remove(2) #He
    #       chemical_numbers.remove(4) #Be
    #       chemical_numbers.remove(45) #Rh
    # QE  chemical_numbers.remove(85) #At
    #     chemical_numbers.remove(18) #Ar
    #     chemical_numbers.remove(53) #I
    # Siesta     chemical_numbers.remove(80) #Hg
    #            chemical_numbers.remove(37) #Rb
    # GPAW     chemical_numbers.remove(83) #Bi
    #          chemical_numbers.remove(43) #Tc
    #          chemical_numbers.remove(86) #Bi
    #          chemical_numbers.remove(84) #Po
    # CP2K    chemical_numbers.remove(11) #Na
    only_must_have_elements = [ase.data.chemical_symbols[i] for i in chemical_numbers]
    
    keep_only_libs = sys.argv[3:]
    if not keep_only_libs:
        keep_only_libs = None

    generate_box_plt('box_plot_cutoffs_',
        only_must_have_elements=only_must_have_elements, material_set_label=material_set_label, file_suffix=elements, keep_only_libs=keep_only_libs)
