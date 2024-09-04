#!/bin/env python

import h5py
import json
from matplotlib import pyplot as plt
import numpy as np
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
    'nc-dojo-v0.4.1-std': '#ffa500',
    'nc-spms-oncvpsp4': '#7f8001', # XXX: new
    "nc-dojo-v0.4.1-str": '#ffb500', #: TBD
    "nc-dojo-v0.5.0-std": '#ffc500', #: TBD
    "nc-sg15-oncvpsp4": '#000000',
    'us-gbrv-v1.x-upf2': '#00cdcd',
    'us-psl-v1.0.0-high': '#ff0000',
    "us-psl-v1.0.0-low": '#fa0000', # TBD
    "us-psl-v0.x": '#0000ff',
    'paw-jth-v1.1-std': '#984ea3', # XXX: new TBD
    "paw-jth-v1.1-str": '#984fa3', # TBD
    "paw-lanthanides-wentzcovitch": '#610b5e',
    "paw-psl-v0.x": '#ff00ff',
    "paw-psl-v1.0.0-high": '#008b00',
    "paw-psl-v1.0.0-low": '#008c00', #TBD
    "paw-actinides-marburg": '#ea388e',
}

lib_abbr_name_mapping = {
    'nc-dojo-v0.4.1-std': 'DOJO-041-std',
    'nc-spms-oncvpsp4': 'SPMS',
    "nc-dojo-v0.4.1-str": 'DOJO-041-str',
    "nc-dojo-v0.5.0-std": 'DOJO-050-std',
    "nc-sg15-oncvpsp4": 'SG15',
    'us-gbrv-v1.x-upf2': 'GBRV-1.X',
    'us-psl-v1.0.0-high': 'PSL-US-v1-high',
    "us-psl-v1.0.0-low": 'PSL-US-v1-low',
    "us-psl-v0.x": 'PSL-US-v0.x',
    'paw-jth-v1.1-std': 'JTH-1.1-std',
    "paw-jth-v1.1-str": 'JTH-1.1-str',
    "paw-lanthanides-wentzcovitch": 'Wentzcovitch',
    "paw-psl-v0.x": 'PSL-PAW-v0.x',
    "paw-psl-v1.0.0-high": 'PSL-PAW-v1-high',
    "paw-psl-v1.0.0-low": 'PSL-PAW-v1-low',
    "paw-actinides-marburg": 'MARBURG',
}

MAX_CUTOFF = 200
MIN_CUTOFF = 30

def plot(element, dual, conff, converge_h5):
    count = 0
    plot_lines = []
    offset = 8

    try:
        pps = element_pps_mapping[element]
    except:
        return

    # Define pyplot instance
    plt.figure(figsize=(40, 4 * len(pps)))

    plt.title(f'Verification summary: {element} ({conff}) (dual={dual}) ({CRITERIA})', fontsize=20)
    plt.xlim(20, 215)
    plt.ylim(-offset/2.,(len(pps)-0.5)*offset)

    # dual = 18 if element in HIGH_DUAL_ELEMENTS else 8
    plt.xlabel(f'Wavefunction cutoff [Ry]; Charge density cutoff [Ry] = {dual} x Ewfc (PAW/US) | 4 x Ewfc (NC); q-point = '+str([0.5, 0.5, 0.5]),fontsize=20)
    plt.ylabel(f'Error w.r.t. ref. wavefunction cutoff (for the SSSP {CRITERIA} criteria)',fontsize=20)

    # legend manually created
    line_phonon_frequencies = Line2D([0], [0], marker='o', linestyle='-', label=r'$\delta \bar{\omega}$', color='black')
    line_cohesive_energy = Line2D([0], [0], marker='v', linestyle='--', label=r'$\delta E_{coh}$', color='black')
    line_pressure = Line2D([0], [0], marker='*', linestyle=':', label=r'$\delta V_{press}$', color='black')
    line_eos = Line2D([0], [0], marker='s', linestyle='-.', label=r'$\delta \nu$', color='black')
    handles = [line_phonon_frequencies, line_cohesive_energy, line_pressure, line_eos]
    legend = plt.legend(handles=handles, bbox_to_anchor=(-0.03, 1.0), fontsize=14, frameon=True, markerscale=1.0)
    plt.gca().add_artist(legend)

    lxticks = list(range(30, 201, 10))
    plt.xticks(lxticks, [str(x) for x in lxticks], fontsize=14)

    # yticks
    ypos = []
    ylabel = []
    for i in range(len(pps)):
        for jp, jl in zip([-2, 0, 2], ['', '0', '']):
            # reset the middle line of every PP to 0
            ypos.append(jp + offset * i)
            ylabel.append(jl)


    plt.yticks(ypos, ylabel)

    for pp_name in pps:
        dataset = converge_h5[pp_name]

        lib_name = dataset.attrs.get('lib_name')
        z_valence = dataset.attrs.get('z_valence')
        if lib_name is None:
            raise ValueError(f"lib_name of {dataset} is None")

        pcolor = lib_color_mapping[lib_name]

        print(f"------> Pseudopotential = {pp_name}")

        # label, Z_valence, TODO: avg.nu, max.nu (conf), friendly-nu (which give lower weight on XO3)
        try:
            eos_dataset = eos_h5[pp_name]
            avg_nu = 0
            n_nu = 0
            avg_nu_wo_xo3 = 0
            n_nu_wo_xo3 = 0
            max_nu = 0
            max_conf = 'n/a'
            for conf, data in eos_dataset['transferability_eos'].items():
                nu = data.attrs.get('nu')
                avg_nu += nu
                n_nu += 1
                if conf != "XO3":
                    avg_nu_wo_xo3 += nu
                    n_nu_wo_xo3 += 1

                # max_nu = max(max_nu, nu)
                if nu > max_nu:
                    max_nu = nu
                    max_conf = conf
            avg_nu /= n_nu
            avg_nu_wo_xo3 /= n_nu_wo_xo3
            text_blob = f"{lib_abbr_name_mapping[lib_name]}\n" + f"Z = {z_valence}\n" + f"avg.$\\nu$ = {avg_nu:.2f}\n" + f"max.$\\nu$ = {max_nu:.2f} ({max_conf})\n" + f"ang.$\\nu$ = {avg_nu_wo_xo3:.2f} (w/o XO3)"
            plt.text(MAX_CUTOFF+21, offset * count, text_blob,
                verticalalignment='center',horizontalalignment='center',fontsize=14)
        except:
            pass

        try:
            xs_phonon_frequencies = dataset['convergence_phonon_frequencies']['xs'][()]
            ys_phonon_frequencies = dataset['convergence_phonon_frequencies']['ys'][()]
            ys_phonon_frequencies *= (2 / PHONON_C_FACTOR)
            ys_phonon_frequencies_max_diff = dataset['convergence_phonon_frequencies']['ys_relative_max_diff'][()]
            low_err_freqs = np.zeros(len(ys_phonon_frequencies))
            high_err_freqs = abs(ys_phonon_frequencies_max_diff) - abs(ys_phonon_frequencies)
            ys_phonon_frequencies += count * offset

            ref_omega_max = dataset['convergence_phonon_frequencies']['ys_omega_max'][-1]
            plt.text(MAX_CUTOFF+13.5, offset*count - 0.6, '$\omega_{max}$ = ' + f'{ref_omega_max:.2f}' + ' cm$^{-1}$',
                horizontalalignment='right',verticalalignment='center',fontsize=14)

            plt.errorbar(xs_phonon_frequencies, ys_phonon_frequencies, yerr=[low_err_freqs, abs(high_err_freqs)], 
                    capthick=3, capsize=4, marker='o', linestyle='-',
                    color=pcolor, alpha=0.8, lw=2, ms=10)        
        except Exception as exc:
            eprint(f"in ploting phonon of {pp_name}: {exc}")

        # Pressure
        try:
            xs_pressure = dataset['convergence_pressure']['xs'][()]
            ys_pressure = dataset['convergence_pressure']['ys'][()]
            ys_pressure *= (2 / PRESSURE_C_FACTOR) # XXX: magnify so the axhline is the criteria, this value is depend on criteria, => 2 / upper_bound(criteria)
            ys_pressure += count * offset
            plt.plot(xs_pressure, ys_pressure, marker='v', linestyle='--',
                        color=pcolor, alpha=0.9, 
                        lw=2, ms=10)        
        except Exception as exc:
            eprint(f"in ploting pressure of {pp_name}: {exc}")

        # Cohesive energy - Heats of formation (e.g. cohesive energies)
        try:
            xs_cohesive_energy = dataset['convergence_cohesive_energy']['xs'][()]
            ys_cohesive_energy = dataset['convergence_cohesive_energy']['ys'][()]
            ys_cohesive_energy += count * offset
            plt.plot(xs_cohesive_energy, ys_cohesive_energy, marker='*', linestyle=':',
                    color=pcolor, alpha=0.9,
                    lw=2,ms=10)
            ref_cohesive_energy_max = dataset['convergence_cohesive_energy']['ys_cohesive_energy_per_atom'][-1]
            plt.text(MAX_CUTOFF+13.5, offset*count + 0.6, '$E_{cov}$ = ' + f'{ref_cohesive_energy_max:.2f}' + ' $meV/atom$',
                horizontalalignment='right',verticalalignment='center',fontsize=14)
        except Exception as exc:
            eprint(f"in ploting cohesive energy of {pp_name}: {exc}")

        # EOS (nu w.r.t 200 Ry)
        try:
            xs_eos = dataset['convergence_eos']['xs'][()]
            ys_eos = dataset['convergence_eos']['ys'][()]
            ys_eos *= 2 / EOS_C_FACTOR # 0.2 is the upper_bound of efficiency criteria
            ys_eos += count * offset
            plt.plot(xs_eos, ys_eos, marker='s', linestyle='-.',
                    color=pcolor, alpha=0.9,
                    lw=2,ms=10)
        except Exception as exc:
            eprint(f"in ploting cohesive energy of {pp_name}: {exc}")

        # Bands data text
        try:
            xs_bands = dataset['convergence_bands']['xs'][()]
            ys_eta_c = dataset['convergence_bands']['ys_eta_c'][()]
            ys_max_diff_c = dataset['convergence_bands']['ys_max_diff_c'][()]

            plt.text(MAX_CUTOFF+5,offset/3+offset*count,"[meV]",
                horizontalalignment='left',verticalalignment='top',fontsize=14)
            plt.text(MAX_CUTOFF+5,-offset/3+offset*count,"[meV]",
                horizontalalignment='left',fontsize=14)
            plt.text(MIN_CUTOFF-2,offset/3+offset*count,'$\eta_{10} =$',horizontalalignment='right',
                verticalalignment='top',fontsize=14)
            plt.text(MIN_CUTOFF-2,-offset/3+offset*count,'$\max \eta_{10} =$',
                horizontalalignment='right',fontsize=14)
            for eta_10, max_10, cutoff in zip(ys_eta_c, ys_max_diff_c, xs_bands):
                if eta_10 > 1000:
                    max_10_text = f"{max_10:.0f}"
                    eta_10_text = f"{eta_10:.0f}"
                else:
                    max_10_text = f"{max_10:.2f}"
                    eta_10_text = f"{eta_10:.2f}"

                plt.text(cutoff, -offset/3+offset*count, max_10_text,
                     color='black',horizontalalignment='center',fontsize=14)
                plt.text(cutoff, offset/3+offset*count, eta_10_text,
                    color='black',horizontalalignment='center',
                    verticalalignment='top',fontsize=14)
        except Exception as exc:
            eprint(f"in ploting bands of {pp_name}: {exc}")
                      
        # increase count so shifting offset for next PP
        count += 1

        # use ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1))) linestyle to seperate
        # https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
        axhlstyle = (0, (3, 1, 1, 1, 1, 1))
        plt.axhline(2+offset*(count-1), 
                color=lib_color_mapping[lib_name], linestyle=axhlstyle, lw=2)  # horizontal line is at y=2
        plt.axhline(-2+offset*(count-1), 
                color=lib_color_mapping[lib_name], linestyle=axhlstyle, lw=2)

        # half precission
        plt.axhline(1+offset*(count-1), 
                color=lib_color_mapping[lib_name], linestyle=(0, (3, 10, 1, 10, 1, 10)), lw=2)  # horizontal line is at y=2


    # plt.savefig(element+'_'+str(dual)+'_conv_patt.png')
    Path(f'plots_{dual}_{CRITERIA}').mkdir(exist_ok=True)
    plt.savefig(f'plots_{dual}_{CRITERIA}/{element}_{conff}_{dual}_{CRITERIA}_summary.png')
    plt.close()

            
if __name__ == "__main__":
    # pre--
    with open('conf_mapping.json', 'r') as fh:
        conf_mapping = json.load(fh)
        conf_mapping = {k: v.lower() for k, v in conf_mapping.items()}

    for conf_dual in ["bcc_dual8", "bcc_dual12", "bcc_dual18", "dc_dual8", "dc_dual12", "dc_dual18"]:
        conff = conf_dual.split('_')[0]
        dual = conf_dual.split('_dual')[1]

        # Load the dataset of convergence results
        converge_h5 = h5py.File(f'./pp_verify_convergence_{conf_dual}.h5')
        eos_h5 = h5py.File('./pp_verify_transferability_eos_200.h5')

        # traverse once to collect mapping of element -> all PPs
        element_pps_mapping = {}

        def curated_by_element(name: str, obj):
            # only get result for first layer
            if '/' in name:
                return
            element = obj.attrs.get('element')
            if element is None:
                raise ValueError(f"element attr of {obj} is None")

            element_pps_mapping.setdefault(element, []).append(name)


        converge_h5.visititems(curated_by_element)

        for element in ALL_ELEMENTS:
            if element == 'Fe':
                import ipdb; ipdb.set_trace()
            # if element and conf not compatible (the stable conf of element), skip
            try:
                stable_conf = conf_mapping[element]
            except:
                continue
            else:
                if stable_conf != conff:
                    continue
            
            print(conff, dual, element)
            plot(element, dual, conff, converge_h5)
