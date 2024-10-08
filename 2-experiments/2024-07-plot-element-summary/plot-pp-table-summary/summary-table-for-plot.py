#!/bin/env python

import h5py
from matplotlib import pyplot as plt
import numpy as np
from collections import defaultdict
from aiida_sssp_workflow.utils.element import ALL_ELEMENTS, HIGH_DUAL_ELEMENTS
import sys
from pathlib import Path
from aiida_sssp_workflow.utils.protocol import get_protocol
import copy
import json

from matplotlib.lines import Line2D

## conf mapping loading
with open('conf_mapping.json', 'r') as fh:
    conf_mapping = json.load(fh)
    conf_mapping = {k: v.lower() for k, v in conf_mapping.items()}

## eos mapping loading
with open('eos.json', 'r') as fh:
    eos_mapping = json.load(fh)

def compute_recommended_cutoffs(xs, ys, criteria) -> int:
    for x, y in zip(reversed(xs), reversed(ys)):
        cutoff = x
        if y > criteria:
            break

    return cutoff

def get_criteria(protocol, property):
    criteria = get_protocol(category='criteria', name=protocol)

    return criteria[property]['bounds'][1]

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Load the dataset of convergence results
eos_h5 = h5py.File('./pp_verify_transferability_eos_200.h5')

converge_h5 = {}
converge_h5['bcc'] = h5py.File(f'./pp_verify_convergence_bcc.h5')
converge_h5['fcc'] = h5py.File(f'./pp_verify_convergence_fcc.h5')
converge_h5['dc'] = h5py.File(f'./pp_verify_convergence_dc.h5')

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


for conf in ['bcc', 'fcc', 'dc']:
    converge_h5['bcc'].visititems(curated_by_element)

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

def compute_rec_cutoff(conff, pp_name, converge_h5): # -> rec_cutoff
    rec_cutoff = {'efficiency': 0, 'precision': 0}
    try:
        dataset = converge_h5[conff][pp_name]
    except:
        rec_cutoff = {'efficiency': None, 'precision': None}
        return rec_cutoff

    lib_name = dataset.attrs.get('lib_name')
    z_valence = dataset.attrs.get('z_valence')
    if lib_name is None:
        raise ValueError(f"lib_name of {dataset} is None")

    # Phonon
    try:
        xs_phonon_frequencies = dataset['convergence_phonon_frequencies']['xs'][()]
        ys_phonon_frequencies = dataset['convergence_phonon_frequencies']['ys'][()]
        cc_ys_phonon_frequencies = copy.copy(ys_phonon_frequencies)

    except Exception as exc:
        eprint(f"in ploting phonon of {pp_name}: {exc}")
        rec_cutoff = {'efficiency': None, 'precision': None}
        return rec_cutoff
    else:
        for protocol in ['efficiency', 'precision']:
            property = "phonon_frequencies"
            criteria = get_criteria(protocol, property)
            co = compute_recommended_cutoffs(xs_phonon_frequencies, cc_ys_phonon_frequencies, criteria)
            if co > rec_cutoff[protocol]:
                rec_cutoff[protocol] = co

    # Pressure
    try:
        xs_pressure = dataset['convergence_pressure']['xs'][()]
        ys_pressure = dataset['convergence_pressure']['ys'][()]
        cc_ys_pressure = copy.copy(ys_pressure)
    except Exception as exc:
        eprint(f"in ploting pressure of {pp_name}: {exc}")
        rec_cutoff = {'efficiency': None, 'precision': None}
        return rec_cutoff
    else:
        for protocol in ['efficiency', 'precision']:
            property = "pressure"
            criteria = get_criteria(protocol, property)
            co = compute_recommended_cutoffs(xs_pressure, cc_ys_pressure, criteria)
            if co > rec_cutoff[protocol]:
                rec_cutoff[protocol] = co

    # Cohesive energy - Heats of formation (e.g. cohesive energies)
    try:
        xs_cohesive_energy = dataset['convergence_cohesive_energy']['xs'][()]
        ys_cohesive_energy = dataset['convergence_cohesive_energy']['ys'][()]
        cc_ys_cohesive_energy = copy.copy(ys_cohesive_energy)
    except Exception as exc:
        eprint(f"in ploting cohesive energy of {pp_name}: {exc}")
        rec_cutoff = {'efficiency': None, 'precision': None}
        return rec_cutoff
    else:
        for protocol in ['efficiency', 'precision']:
            property = "cohesive_energy"
            criteria = get_criteria(protocol, property)
            co = compute_recommended_cutoffs(xs_cohesive_energy, cc_ys_cohesive_energy, criteria)
            if co > rec_cutoff[protocol]:
                rec_cutoff[protocol] = co


    # EOS (nu w.r.t 200 Ry)
    try:
        xs_eos = dataset['convergence_eos']['xs'][()]
        ys_eos = dataset['convergence_eos']['ys'][()]
        cc_ys_eos = copy.copy(ys_eos)
    except Exception as exc:
        eprint(f"in ploting cohesive energy of {pp_name}: {exc}")
        rec_cutoff = {'efficiency': None, 'precision': None}
        return rec_cutoff
    else:
        for protocol in ['efficiency', 'precision']:
            property = "eos"
            criteria = get_criteria(protocol, property)
            co = compute_recommended_cutoffs(xs_eos, cc_ys_eos, criteria)
            if co > rec_cutoff[protocol]:
                rec_cutoff[protocol] = co

    # Bands data text
    try:
        xs_bands = dataset['convergence_bands']['xs'][()]
        ys_eta_c = dataset['convergence_bands']['ys_eta_c'][()]
        cc_ys_eta_c = copy.copy(ys_eta_c)
    except Exception as exc:
        eprint(f"in ploting bands of {pp_name}: {exc}")
        rec_cutoff = {'efficiency': None, 'precision': None}
        return rec_cutoff
    else:
        for protocol in ['efficiency', 'precision']:
            property = "bands"
            criteria = get_criteria(protocol, property)
            co = compute_recommended_cutoffs(xs_bands, cc_ys_eta_c, criteria)
            if co > rec_cutoff[protocol]:
                rec_cutoff[protocol] = co

    return rec_cutoff

def summary(element):
    ele_full_report = {}

    try:
        pps = element_pps_mapping[element]
    except:
        return {}

    stable_conf = conf_mapping[element]

    for pp_name in pps:
        ereport = {}
        ## for test 
        #if not 'us' in pp_name:
        #    continue
        print(f"------> Pseudopotential = {pp_name}")

        ereport['element'] = element
        ereport['stable_conf'] = stable_conf

        z_valence = pp_name.split(".")[3].split("_")[1]
        ereport['z_valence'] = z_valence

        # label, Z_valence, TODO: avg.nu, max.nu (conf), friendly-nu (which give lower weight on XO3)
        try:
            eos_dataset = eos_h5[pp_name]
            tot_nu = 0
            n_nu = 0
            avg_nu_wo_max = 0
            n_nu_wo_max = 0
            max_nu = 0
            max_conf = 'n/a'

            nu_confs = {}
            for conf, data in eos_dataset['transferability_eos'].items():
                nu = data.attrs.get('nu')
                tot_nu += nu
                n_nu += 1

                nu_confs[conf] = nu

                # max_nu = max(max_nu, nu)
                if nu > max_nu:
                    max_nu = nu
                    max_conf = conf

            avg_nu = tot_nu / n_nu
            avg_nu_wo_max = (tot_nu - max_nu) / (n_nu - 1)

            ereport['nu_confs'] = nu_confs

            # ereport['avg.nu'] = avg_nu
            ereport['avg.nu.w/o.max'] = avg_nu_wo_max
            ereport['nu_eff_score'] = compute_nu_eff_score(eos_dataset['transferability_eos'])
            ereport['nu_prec_score'] = compute_nu_prec_score(eos_dataset['transferability_eos'])
        except:
            #import ipdb; ipdb.set_trace()
            # ereport['avg.nu'] = None
            ereport['avg.nu.w/o.max'] = None
            ereport['nu_eff_score'] = None
            ereport['nu_prec_score'] = None

        rec_cutoff = compute_rec_cutoff(stable_conf, pp_name, converge_h5)      

        for cri in ['efficiency', 'precision']:
            ereport[cri] = rec_cutoff[cri]

            # end of conff loop

        ele_full_report[pp_name] = ereport

    # end
    return ele_full_report

def compute_nu_eff_score(nu_data):
    base_nu = 0.33
    
    max_nu = 0
    tot_score = 0
    n_nu = 0
    for conf, data in nu_data.items():
        nu = data.attrs.get('nu')
        if nu > base_nu:
            tot_score += nu - base_nu

        n_nu += 1

		# max_nu = max(max_nu, nu)
        if nu > max_nu:
            max_nu = nu
            max_conf = conf
    tot_score = tot_score - (max_nu - base_nu)
    score = tot_score / (n_nu - 1)
    return score

def compute_nu_prec_score(nu_data):
    base_nu = 0.1
    
    max_nu = 0
    tot_score = 0
    n_nu = 0
    for conf, data in nu_data.items():
        nu = data.attrs.get('nu')
        if nu > base_nu:
            tot_score += nu - base_nu

        n_nu += 1

		# max_nu = max(max_nu, nu)
        if nu > max_nu:
            max_nu = nu
            max_conf = conf
    tot_score = tot_score - (max_nu - base_nu)
    score = tot_score / (n_nu - 1)
    return score

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
            
if __name__ == "__main__":

    all_report = {}
    # for element in ALL_ELEMENTS:
    for element in ['H', 'Fe']:
        report = summary(element)
        all_report.update(report)

    with open('table-report-tmp.json', 'w') as outfile:
        json.dump(all_report, fp=outfile, cls=NpEncoder, indent=4)
