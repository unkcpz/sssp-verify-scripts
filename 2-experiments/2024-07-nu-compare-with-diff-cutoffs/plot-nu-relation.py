#!/bin/env python
import matplotlib.pyplot as plt
import sys
import numpy as np
import h5py

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

f200 = h5py.File('./pp_verify_transferability_eos_200.h5')
fprec = h5py.File('./pp_verify_transferability_eos_precision.h5')
feff = h5py.File('./pp_verify_transferability_eos_efficiency.h5')

# traverse once to collect mapping of element -> all PPs
lib_pps_mapping = {}

def curated_by_lib(name: str, obj):
    # only get result for first layer
    if '/' in name:
        return
    lib_name = obj.attrs.get('lib_name')
    if lib_name is None:
        raise ValueError(f"lib_name attr of {obj} is None")

    lib_pps_mapping.setdefault(lib_name, []).append(name)


f200.visititems(curated_by_lib)

# get data
x_eff = []
y_eff = []
x_prec = []
y_prec = []
# lib_name = 'nc-dojo-v0.4.1-std'
for lib_name in [
    'nc-dojo-v0.4.1-std',
    'paw-jth-v1.1-std',
    'us-gbrv-v1.x-upf2',
]:
    pps = lib_pps_mapping[lib_name]

    for pp_name in pps:
        try:
            eos_dataset_200 = f200[pp_name]
            eos_dataset_eff = feff[pp_name]
            for conf in eos_dataset_200['transferability_eos'].keys():
                data_200 = eos_dataset_200['transferability_eos'][conf]
                nu_200 = data_200.attrs.get('nu')

                data_eff = eos_dataset_eff['transferability_eos'][conf]
                nu_eff = data_eff.attrs.get('nu')

                x_eff.append(nu_200)
                y_eff.append(nu_eff)
                if abs(nu_eff - nu_200) > 0.05:
                    print(f"eff: {conf} of {pp_name}, nu: {nu_eff:.2f}, nu_200: {nu_200:.2f}")
        except KeyError as kexc:
            eprint(f'{kexc}')

        try:
            eos_dataset_200 = f200[pp_name]
            eos_dataset_prec = fprec[pp_name]
            for conf in eos_dataset_200['transferability_eos'].keys():
                data_200 = eos_dataset_200['transferability_eos'][conf]
                nu_200 = data_200.attrs.get('nu')

                data_prec = eos_dataset_prec['transferability_eos'][conf]
                nu_prec = data_prec.attrs.get('nu')

                x_prec.append(nu_200)
                y_prec.append(nu_prec)
                if abs(nu_prec - nu_200) > 0.05:
                    print(f"prec: {conf} of {pp_name}, nu: {nu_prec:.2f}, nu_200: {nu_200:.2f}")
        except KeyError as kexc:
            eprint(f'{kexc}')

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
fig.subplots_adjust(wspace=0.3)
fig.suptitle(r'Relation between $\nu$ from 200 Ry / from rec cutoffs')
axs[0].scatter(x_eff, y_eff, marker='s', color='blue', label='Effeciecy')
axs[0].scatter(x_prec, y_prec, marker='P', color='orange', label='precision')
axs[0].legend()
axs[0].set_xlim([-0.1, 0.5])
axs[0].set_ylim([-0.1, 0.5])
axs[0].set_xlabel(r"$\nu$ at (200, 200 * dual) Ry cutoffs")
axs[0].set_ylabel(r"$\nu$ at 'converged' cutoffs")

# yyd_eff = (np.array(y_eff) - np.array(x_eff)) / (np.array(x_eff))
yyd_eff = np.array(y_eff) - np.array(x_eff)
xxd_eff = np.arange(start=0, stop=len(yyd_eff))
# yyd_prec = (np.array(y_prec) - np.array(x_prec)) / (np.array(x_prec))
yyd_prec = np.array(y_prec) - np.array(x_prec)
xxd_prec = np.arange(start=0, stop=len(yyd_prec))
axs[1].scatter(xxd_eff, yyd_eff, marker='s', color='blue', label='Effeciency')
axs[1].scatter(xxd_prec, yyd_prec, marker='P', color='orange', label='precision')
axs[1].set_xticks([])
axs[1].set_xlabel("All verifications")
axs[1].set_ylabel(r"abs($\nu_{200ry}$ - $\nu_{c}$)")

# yyd_eff_ex = []
# for x, y in zip(x_eff, y_eff):
#     if x > 1.0:
#         continue
#     yyd_eff_ex.append(y - x)
# xxd_eff_ex = np.arange(start=0, stop=len(yyd_eff_ex))
#
# yyd_prec_ex = []
# for x, y in zip(x_prec, y_prec):
#     if x > 1.0:
#         continue
#     yyd_prec_ex.append(y - x)
# xxd_prec_ex = np.arange(start=0, stop=len(yyd_prec_ex))
# axs[2].scatter(xxd_eff_ex, yyd_eff_ex, marker='s', color='blue', label='Effeciency')
# axs[2].scatter(xxd_prec_ex, yyd_prec_ex, marker='P', color='orange', label='precision')

plt.savefig('eos-relation.png')
plt.close()
