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
}

pk_name_mapping = {
    "4874271": "Ni.nc.pbe.z_18.oncvpsp3.dojo.v0.4.1-str.upf",
    "4874280": "Ni.us.pbe.z_18.uspp.gbrv.v1.4.upf",
    "4874294": "Ni.us.pbe.z_10.ld1.psl.v1.0.0-low.upf",
    "4874312": "Ni.paw.pbe.z_18.ld1.psl.v1.0.0-high.upf",
    "4874333": "Ni.us.pbe.z_10.ld1.psl.v0.1.upf",
    "4874373": "Ni.nc.pbe.z_18.oncvpsp4.spms.v1.upf",
    "4941897": "Ni.paw.pbe.z_10.ld1.psl.v1.0.0-low.upf",
    "4941906": "Ni.nc.pbe.z_18.oncvpsp4.sg15.v0.upf",
    "4941919": "Ni.paw.pbe.z_10.ld1.psl.v0.1.upf",
}

# 4955804  Co.nc.pbe.z_17.oncvpsp4.sg15.v0.upf                      0
# 4955814  Co.us.pbe.z_17.ld1.psl.v1.0.0-high.upf                   0
# 4955829  Co.nc.pbe.z_17.oncvpsp3.dojo.v0.4.1-std.upf              0
# 4955844  Co.paw.pbe.z_17.atompaw.jth.v1.1-std.upf                 0
# 4955868  Co.us.pbe.z_17.uspp.gbrv.v1.2.upf                        0
# 4955907  Co.paw.pbe.z_9.ld1.psl.v1.0.0-low.upf                    0
# 4983614  Co.nc.pbe.z_17.oncvpsp4.spms.v1.upf                      0
# 4983624  Co.us.pbe.z_9.ld1.psl.v1.0.0-low.upf                     0
# 4983637  Co.paw.pbe.z_17.ld1.psl.v1.0.0-high.upf                  0
# 4983655  Co.nc.pbe.z_17.oncvpsp3.dojo.v0.4.1-str.upf              0


def find_lib(name):
    if "us" in name and "psl.v0" in name:
        return "us-psl-v0.x"
    elif "us" in name and "psl.v1.0.0-high" in name:
        return "us-psl-v1.0.0-high"
    elif "us" in name and "psl.v1.0.0-low" in name:
        return "us-psl-v1.0.0-low"
    elif "us" in name and "gbrv" in name:
        return "us-gbrv-v1.x-upf2"
    elif "paw" in name and "psl.v0" in name:
        return "paw-psl-v0.x"
    elif "paw" in name and "psl.v1.0.0-high" in name:
        return "paw-psl-v1.0.0-high"
    elif "paw" in name and "psl.v1.0.0-low" in name:
        return "paw-psl-v1.0.0-low"
    elif "paw" in name and "jth.v1.1-std" in name:
        return "paw-jth-v1.1-std"
    elif "paw" in name and "jth.v1.1-str" in name:
        return "paw-jth-v1.1-str"
    elif "nc" in name and "sg15" in name:
        return "nc-sg15-oncvpsp4"
    elif "nc" in name and "spms" in name:
        return "nc-spms-oncvpsp4"
    elif "nc" in name and "dojo.v0.4.1-std" in name:
        return "nc-dojo-v0.4.1-std"
    elif "nc" in name and "dojo.v0.4.1-str" in name:
        return "nc-dojo-v0.4.1-str"
    else:
        print(name)
        return "unknown pp"


# %%
fig, ax = plt.subplots()
balence_v = 10.99

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
            x = 100 * (volume - balence_v) / balence_v

            if volume > 50:
                break

            xs.append(x)
            ys.append(mag)
        except:
            pass
    lib_name = find_lib(name)
    label = lib_abbr_name_mapping[lib_name]
    color = lib_color_mapping[lib_name]
    ax.plot(xs, ys, label=label, color=color, marker="o")


# %%
ax.set_xlim(-20, 10)
ax.set_ylim(0, 0.8)
ax.set_xlabel("Volume variation (%)")
ax.set_ylabel("magnetization (muB)")
ax.legend(loc=2, prop={"size": 8})
fig.set_size_inches((10, 8))
fig.savefig("/tmp/test-ni.png")

# %%

