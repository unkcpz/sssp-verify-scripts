#!/bin/env python

import h5py
from matplotlib import pyplot as plt
import numpy as np
import json
from aiida_sssp_workflow.utils.element import ALL_ELEMENTS, HIGH_DUAL_ELEMENTS
import sys
from pathlib import Path

from matplotlib.lines import Line2D


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# cri
CRITERIA = "efficiency"

if CRITERIA == "efficiency":
    EOS_C_FACTOR = 0.2
    PRESSURE_C_FACTOR = 1
    PHONON_C_FACTOR = 2
elif CRITERIA == "precision":
    EOS_C_FACTOR = 0.1
    PRESSURE_C_FACTOR = 0.5
    PHONON_C_FACTOR = 1

# pseudos_colors_dict = dict([(pseudo,color) for pseudo,color in zip(
#              ['100PAW','100PAW_low','100US','100US_low','031PAW','031US',
#               'GBRV-1.2','GBRV-1.4','GBRV-1.5','SG15','SG15-1.1','Goedecker',
#               'Dojo','THEOS','Wentzcovitch','Vanderbilt','THEOS2','all_elec',
#               'all_elec_denser','BM','GIPAW','psorigPAW','psorigUS'],
#              ['#008B00','#80FF80','#FF0000','#FF8080','#FF00FF','#0000FF',
#               '#00CDCD','#4682B4','#B8860B','#000000','#708090','#808000',
#               '#FFA500','#D7DF01','#610B5E','#8FBC8F','#F0F000','#F000F0',
#               '#00F0F0','#A5FF00','#B44682','#CD00CD','#86B80B']
#               )])
lib_color_mapping = {
    "nc-dojo-v0.4.1-std": "#ffa500",
    "nc-spms-oncvpsp4": "#7f8001",  # XXX: new
    "nc-dojo-v0.4.1-str": "#ffb500",  #: TBD
    "nc-dojo-v0.5.0-std": "#ffc500",  #: TBD
    "nc-sg15-oncvpsp4": "#000000",
    "us-gbrv-v1.x-upf2": "#00cdcd",
    "us-psl-v1.0.0-high": "#ff0000",
    "us-psl-v1.0.0-low": "#fa0000",  # TBD
    "us-psl-v0.x": "#0000ff",
    "paw-jth-v1.1-std": "#984ea3",  # XXX: new TBD
    "paw-jth-v1.1-str": "#984fa3",  # TBD
    "paw-lanthanides-wentzcovitch": "#610b5e",
    "paw-psl-v0.x": "#ff00ff",
    "paw-psl-v1.0.0-high": "#008b00",
    "paw-psl-v1.0.0-low": "#008c00",  # TBD
    "paw-actinides-marburg": "#ea388e",
    "sssp-eff-v2-0911": "#ffa500",
}

lib_abbr_name_mapping = {
    "nc-dojo-v0.4.1-std": "DOJO-041-std",
    "nc-spms-oncvpsp4": "SPMS",
    "nc-dojo-v0.4.1-str": "DOJO-041-str",
    "nc-dojo-v0.5.0-std": "DOJO-050-std",
    "nc-sg15-oncvpsp4": "SG15",
    "us-gbrv-v1.x-upf2": "GBRV-1.X",
    "us-psl-v1.0.0-high": "PSL-US-v1-high",
    "us-psl-v1.0.0-low": "PSL-US-v1-low",
    "us-psl-v0.x": "PSL-US-v0.x",
    "paw-jth-v1.1-std": "JTH-1.1-std",
    "paw-jth-v1.1-str": "JTH-1.1-str",
    "paw-lanthanides-wentzcovitch": "Wentzcovitch",
    "paw-psl-v0.x": "PSL-PAW-v0.x",
    "paw-psl-v1.0.0-high": "PSL-PAW-v1-high",
    "paw-psl-v1.0.0-low": "PSL-PAW-v1-low",
    "paw-actinides-marburg": "MARBURG",
    "sssp-eff-v2-0911": "SSSP-v2-eff",
}

MAX_CUTOFF = 200
MIN_CUTOFF = 30

prop_color_mapping = {
    "phonon": "#008b00",
    "cohesive": "#ea388e",
    "pressure": "#0000ff",
    "eos": "#984ea3",
    "bands_eta": "#00cdcd",
    "bands_max": "#00000a",
}


def plot(ax, element, conff, converge_h5):
    try:
        pps = element_pps_mapping[element]
    except:
        # raise
        return

    # --------------------------------------------------------
    pp_name = pps[0]
    dataset = converge_h5[pp_name]

    # Phonon
    try:
        xs_phonon = dataset["convergence_phonon_frequencies"]["xs"][()]
        ys_phonon = dataset["convergence_phonon_frequencies"]["ys"][()]
        ys_phonon *= 2 / PHONON_C_FACTOR
        ax.plot(
            xs_phonon,
            ys_phonon,
            marker="o",
            linestyle="-",
            color=prop_color_mapping["phonon"],
        )
    except Exception as exc:
        eprint(f"in plotting phonon of {pp_name}: {exc}")

    # Pressure
    try:
        xs_pressure = dataset["convergence_pressure"]["xs"][()]
        ys_pressure = dataset["convergence_pressure"]["ys"][()]
        ys_pressure *= 2 / PRESSURE_C_FACTOR
        ax.plot(
            xs_pressure,
            ys_pressure,
            marker="v",
            linestyle="--",
            color=prop_color_mapping["pressure"],
            alpha=0.9,
            lw=2,
            ms=10,
        )
    except Exception as exc:
        eprint(f"in plotting pressure of {pp_name}: {exc}")

    # Cohesive energy
    try:
        xs_cohesive_energy = dataset["convergence_cohesive_energy"]["xs"][()]
        ys_cohesive_energy = dataset["convergence_cohesive_energy"]["ys"][()]
        ax.plot(
            xs_cohesive_energy,
            ys_cohesive_energy,
            marker="*",
            linestyle=":",
            color=prop_color_mapping["cohesive"],
            alpha=0.9,
            lw=2,
            ms=10,
        )
        # ref_cohesive_energy_max = dataset['convergence_cohesive_energy']['ys_cohesive_energy_per_atom'][-1]
        # ax.text(MAX_CUTOFF + 13.5,  0.6,
        #         '$E_{cov}$ = ' + f'{ref_cohesive_energy_max:.2f}' + ' $meV/atom$',
        #         ha='right', va='center', fontsize=14)
    except Exception as exc:
        ax.text(
            MAX_CUTOFF + 13.5,
            0.6,
            "$E_{cov}$ (not avail) ",
            ha="right",
            va="center",
            fontsize=14,
        )
        eprint(f"in plotting cohesive energy of {pp_name}: {exc}")

    # EOS
    try:
        xs_eos = dataset["convergence_eos"]["xs"][()]
        ys_eos = dataset["convergence_eos"]["ys"][()]
        ys_eos *= 2 / EOS_C_FACTOR
        ax.plot(
            xs_eos,
            ys_eos,
            marker="s",
            linestyle="-.",
            color=prop_color_mapping["eos"],
            alpha=0.9,
            lw=2,
            ms=10,
        )
    except Exception as exc:
        eprint(f"in plotting EOS of {pp_name}: {exc}")

    # Bands
    try:
        xs_bands = dataset["convergence_bands"]["xs"][()]
        ys_eta_c = dataset["convergence_bands"]["ys_eta_c"][()]
        ys_max_diff_c = dataset["convergence_bands"]["ys_max_diff_c"][()]

        # ax.text(MAX_CUTOFF+5, 0, "[meV]",
        #         ha='left', va='top', fontsize=14)
        # ax.text(MAX_CUTOFF+5, 0, "[meV]",
        #         ha='left', fontsize=14)
        # ax.text(MIN_CUTOFF-2, 0, '$\eta_{10} =$',
        #         ha='right', va='top', fontsize=14)
        # ax.text(MIN_CUTOFF-2, 0, '$\max \eta_{10} =$',
        #         ha='right', fontsize=14)
        #
        for eta_10, max_10, cutoff in zip(ys_eta_c, ys_max_diff_c, xs_bands):
            if eta_10 > 1000:
                max_10_text = f"{max_10:.0f}"
                eta_10_text = f"{eta_10:.0f}"
            else:
                max_10_text = f"{max_10:.2f}"
                eta_10_text = f"{eta_10:.2f}"

            ax.text(
                cutoff,
                -0.35,
                eta_10_text,
                color=prop_color_mapping["bands_eta"],
                ha="center",
                fontsize=14,
            )
            ax.text(
                cutoff,
                -0.15,
                max_10_text,
                color=prop_color_mapping["bands_max"],
                ha="center",
                fontsize=14,
            )
    except Exception as exc:
        eprint(f"in plotting bands of {pp_name}: {exc}")
    # --------------------------------------------------------

    try:
        xs = dataset["convergence_phonon_frequencies"]["xs"][()]
    except:
        # xs = dataset["convergence_eos"]["xs"][()]
        return


    xs_max = max(xs)

    pp_name: str
    pp_type = pp_name.split(".")[1]
    if pp_type == "nc":
        dual = 4

    if pp_type == "paw" or pp_type == "us":
        dual = 8

    if (pp_type == "paw" or pp_type == "us") and element in HIGH_DUAL_ELEMENTS:
        dual = 18

    ecutwfc = int(xs_max / dual)
    ax.set_xlabel(
        f"Charge density cutoff [Ry];  Wavefunction cutoff [Ry] = {ecutwfc}",
        fontsize=12,
    )
    ax.set_ylabel(f"Error w.r.t. ref. cutoff ({CRITERIA})", fontsize=12)

    # yticks setup
    ypos, ylabel = [], []
    for jp, jl in zip([-2, 0, 2], ["", "0", ""]):
        ypos.append(jp)
        ylabel.append(jl)
    ax.set_yticks(ypos)
    ax.set_yticklabels(ylabel)

    ax.set_xticks(xs)

    # make sure limits reflect the plotted artists
    ax.relim()
    ax.autoscale_view()

    ymin, ymax = ax.get_ylim()
    if ymin < -10 or ymax > 10:
        ax.set_ylim(-10, 10)

    ax.set_yticks(ypos)
    ax.set_yticklabels(ylabel)
    ax.set_xticks(xs)


    axhlstyle = (0, (3, 1, 1, 1, 1, 1))

    ax.axhline(2, color="black", linestyle=axhlstyle, lw=2)

    ax.axhline(-2, color="black", linestyle=axhlstyle, lw=2)

    # half precision
    ax.axhline(
        1,
        color="black",
        linestyle=(0, (3, 10, 1, 10, 1, 10)),
        lw=2,
    )

    ax.set_title(f"{pp_name} ({conff}) ({CRITERIA})", fontsize=14)

    return


if __name__ == "__main__":
    from tqdm import tqdm
    import os

    os.makedirs("rho_conv_plots", exist_ok=True)

    with open("conf_mapping.json", "r") as fh:
        conf_mapping = json.load(fh)
        conf_mapping = {k: v.lower() for k, v in conf_mapping.items()}

    pairs = []

    converge_h5 = h5py.File(f"./pp_verify_convergence_final_conf.h5")

    selected_elements = [el for el in ALL_ELEMENTS]

    for el in selected_elements:
        pairs.append((el, converge_h5))

    element_pps_mapping = {}

    def curated_by_element(name: str, obj):
        if "/" in name:
            return
        element = obj.attrs.get("element")
        if element is None:
            raise ValueError(f"element attr of {obj} is None")

        element_pps_mapping.setdefault(element, []).append(name)

    converge_h5.visititems(curated_by_element)

    # now do one figure per element
    for element, converge_h5 in tqdm(pairs, total=len(pairs)):
        fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=False)

        # call your plot() on this one axis
        with open("./conf_mapping.json", mode="r") as fh:
            conf_map: dict[str, str] = json.load(fh)

        try:
            conf = conf_map[element]
        except:
            continue

        plot(ax, element, conf, converge_h5)

        # shared legend on each figure
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="-",
                label=r"$\delta \bar{\omega}$",
                color=prop_color_mapping["phonon"],
            ),
            Line2D(
                [0],
                [0],
                marker="v",
                linestyle="--",
                label=r"$\delta E_{coh}$",
                color=prop_color_mapping["cohesive"],
            ),
            Line2D(
                [0],
                [0],
                marker="*",
                linestyle=":",
                label=r"$\delta V_{press}$",
                color=prop_color_mapping["pressure"],
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="-.",
                label=r"$\delta \nu$",
                color=prop_color_mapping["eos"],
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="-.",
                label=r"$\delta \eta$",
                color=prop_color_mapping["bands_eta"],
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="-.",
                label=r"$\delta \eta_{max}$",
                color=prop_color_mapping["bands_max"],
            ),
        ]
        fig.legend(
            handles=handles, loc="upper center", ncol=len(handles), fontsize=12, frameon=True
        )

        # save into folder
        out_path = os.path.join("rho_conv_plots", f"{element}.png")
        fig.savefig(out_path)
        plt.close(fig)
