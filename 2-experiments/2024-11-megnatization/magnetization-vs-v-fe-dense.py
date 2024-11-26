# %%
from aiida import orm
import aiida
import matplotlib.pyplot as plt

aiida.load_profile()

# %%
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
    "us-gbrv-v1.x-upf2-FD": "#00cdcd",
    "us-psl-v0.x-FD": "#0000ff",
    "us-gbrv-v1.x-upf2-COLD": "#00cdcd",
    "us-psl-v0.x-COLD": "#0000ff",
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
    "us-gbrv-v1.x-upf2-FD": "GBRV-1.X-FD",
    "us-psl-v0.x-FD": "PSL-US-v0.x-FD",
    "us-gbrv-v1.x-upf2-COLD": "GBRV-1.X-COLD",
    "us-psl-v0.x-COLD": "PSL-US-v0.x-COLD",
}

# pk_name_mapping = {
#     "4786669": "Fe.us.pbe.z_16.ld1.psl.v0.2.1.upf",
#     "4786678": "Fe.us.pbe.z_16.ld1.psl.v1.0.0-high.upf",
#     "4786691": "Fe.nc.pbe.z_16.oncvpsp4.sg15.v0.upf",
#     "4786703": "Fe.paw.pbe.z_8.ld1.psl.v1.0.0-low.upf",
#     "4806038": "Fe.paw.pbe.z_16.ld1.psl.v1.0.0-high.upf",
#     "4806047": "Fe.paw.pbe.z_16.atompaw.jth.v1.1-std.upf",
#     "4806059": "Fe.us.pbe.z_8.ld1.psl.v1.0.0-low.upf",
#     "4806075": "Fe.us.pbe.z_16.uspp.gbrv.v1.5.upf",
#     "4824369": "Fe.paw.pbe.z_16.ld1.psl.v0.2.1.upf",
#     "4824379": "Fe.nc.pbe.z_16.oncvpsp4.spms.v1.upf",
#     "4824392": "Fe.nc.pbe.z_16.oncvpsp3.dojo.v0.4.1-std.upf",
#     "4824414": "Fe.nc.pbe.z_16.oncvpsp3.dojo.v0.4.1-str.upf",
# }

pk_name_mapping = {
    "4806075": "Fe.us.pbe.z_16.uspp.gbrv.v1.5-fd.upf",
    "4786669": "Fe.us.pbe.z_16.ld1.psl.v0.2.1-fd.upf",
    "4932627": "Fe.us.pbe.z_16.uspp.gbrv.v1.5-cold.upf",
    "4932636": "Fe.us.pbe.z_16.ld1.psl.v0.2.1-cold.upf",
}


def find_lib(name):
    # if "us" in name and "psl.v0.2.1" in name:
    #     return "us-psl-v0.x"
    # elif "us" in name and "psl.v1.0.0-high" in name:
    #     return "us-psl-v1.0.0-high"
    # elif "us" in name and "psl.v1.0.0-low" in name:
    #     return "us-psl-v1.0.0-low"
    # elif "us" in name and "gbrv" in name:
    #     return "us-gbrv-v1.x-upf2"
    # elif "paw" in name and "psl.v0.2.1" in name:
    #     return "paw-psl-v0.x"
    # elif "paw" in name and "psl.v1.0.0-high" in name:
    #     return "paw-psl-v1.0.0-high"
    # elif "paw" in name and "psl.v1.0.0-low" in name:
    #     return "paw-psl-v1.0.0-low"
    # elif "paw" in name and "jth.v1.1-std" in name:
    #     return "paw-jth-v1.1-std"
    # elif "paw" in name and "jth.v1.1-str" in name:
    #     return "paw-jth-v1.1-str"
    # elif "nc" in name and "sg15" in name:
    #     return "nc-sg15-oncvpsp4"
    # elif "nc" in name and "spms" in name:
    #     return "nc-spms-oncvpsp4"
    # elif "nc" in name and "dojo.v0.4.1-std" in name:
    #     return "nc-dojo-v0.4.1-std"
    # elif "nc" in name and "dojo.v0.4.1-str" in name:
    #     return "nc-dojo-v0.4.1-str"
    # else:
    #     print(name)
    #     return "unknown pp"
    if "us" in name and "psl.v0.2.1" in name and "fd" in name:
        return "us-psl-v0.x-FD"
    elif "us" in name and "gbrv" in name and "fd" in name:
        return "us-gbrv-v1.x-upf2-FD"
    elif "us" in name and "psl.v0.2.1" in name and "cold" in name:
        return "us-psl-v0.x-COLD"
    elif "us" in name and "gbrv" in name and "cold" in name:
        return "us-gbrv-v1.x-upf2-COLD"
    else:
        print(name)
        return "unknown pp"


# %%
fig, ax = plt.subplots()
fe_v = 11.38

# %%
for pk, name in pk_name_mapping.items():
    node = orm.load_node(pk)
    xs = []
    ys = []
    for report in node.outputs.report.get_dict()["report_list"]:
        uuid = report["uuid"]
        n = orm.load_node(uuid)

        try:
            structure = n.inputs.structure
            mag = n.outputs.output_parameters.get_dict()["total_magnetization"]

            volume = structure.get_cell_volume()
            x = 100 * (volume - fe_v) / fe_v

            if volume > 50:
                break

            xs.append(x)
            ys.append(mag)
        except:
            pass

    lib_name = find_lib(name)
    label = lib_abbr_name_mapping[lib_name]
    color = lib_color_mapping[lib_name]
    if "FD" in label:
        marker = "o"
    else:
        marker = "x"
    ax.plot(xs, ys, label=label, color=color, marker=marker)

# %%
ax.set_xlim(-5, 12)
ax.set_ylim(2.1, 2.7)
ax.set_xlabel("Volume variation (%)")
ax.set_ylabel("magnetization (muB)")
ax.legend(loc=2, prop={"size": 8})
fig.set_size_inches((10, 8))
fig.savefig("/tmp/test.png")

# %%
