import json

with open("./cutoffs-orig.json", mode="r") as f:
    map = json.load(f)
    map = {k: v["cutoff_wfc"] for k, v in map.items()}

with open("cutoff-only-ecutwfc.json", mode="w") as f:
    json.dump(map, f, indent=2)


