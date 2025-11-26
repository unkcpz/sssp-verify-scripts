import json


def run():
    with open("./cutoffs-orig.json", mode="r") as fh:
        odata = json.load(fh)

    with open("./rho_conv.json", mode="r") as fh:
        rhodata = json.load(fh)

    new_data = {}
    for k, v in odata.items():
        new_rho = rhodata.get(k)
        if new_rho:
            new_data[k] = {"cutoff_wfc": v["cutoff_wfc"], "cutoff_rho": new_rho}
        else:
            new_data[k] = v

    with open("./cutoff_rho_converge.json", mode="w") as fh:
        json.dump(new_data, fh, indent=4)

if __name__ == "__main__":
    run()
