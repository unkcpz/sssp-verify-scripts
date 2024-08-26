import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.patches import Patch

# The dictionary with the data
with open('table-report.json', 'r') as fh:
	data = json.load(fh)

# Extracting data
rows = []
for key, value in data.items():
    ffn = '.'.join(key.split('.')[:-1])
    rows.append([ffn, 
        round(value["avg.nu"], 2) if value["avg.nu"] is not None else "N/A",
        round(value["avg.nu.w/o.max"], 2) if value["avg.nu.w/o.max"] is not None else "N/A",
        value["bcc.efficiency"] if value["bcc.efficiency"] is not None else "N/A",
        value["bcc.precision"] if value["bcc.precision"] is not None else "N/A",
        value["fcc.efficiency"] if value["fcc.efficiency"] is not None else "N/A",
        value["fcc.precision"] if value["fcc.precision"] is not None else "N/A",
        value["dc.efficiency"] if value["dc.efficiency"] is not None else "N/A",
        value["dc.precision"] if value["dc.precision"] is not None else "N/A",
        round(value["eff_score"], 2) if value["eff_score"] is not None else "N/A",
        round(value["prec_score"], 2) if value["prec_score"] is not None else "N/A"],
     )

# Convert to numpy array for easier handling
rows = np.array(rows, dtype=object)

# Group by element and insert blank rows
element_previous = rows[0, 0].split('.')[0]  # Initial element
rows_with_blanks = []
min_indices = []
for i, row in enumerate(rows):
    element_current = row[0].split('.')[0]  # Current element
    if element_current != element_previous:
        # Insert a blank row when the element changes
        rows_with_blanks.append([""] * len(row))
        element_previous = element_current
    rows_with_blanks.append(row)

# Convert back to a numpy array
rows = np.array(rows_with_blanks, dtype=object)

# Group by element
elements = {}
for i, row in enumerate(rows):
    element = row[0].split('.')[0]
    if element not in elements:
        elements[element] = []
    elements[element].append(i)

# Create the table plot
fig, ax = plt.subplots(figsize=(30, 8))
ax.axis('tight')
ax.axis('off')

# Create the table
# Define column widths
col_widths = [0.15] + [0.05] * (len(rows[0]) - 1)
table = ax.table(cellText=rows, colLabels=["Pseudo", "avg.nu", "avg.nu w/o max",
                                           "bcc eff (Ry)", "bcc prec (Ry)",
                                           "fcc eff (Ry)", "fcc prec (Ry)",
                                           "dc eff (Ry)", "dc prec (Ry)", 
                                           "Eff score", "Prec score"],
                 cellLoc='center', loc='center', colWidths=col_widths)

# Find the index of the row with the minimum avg.nu for each element
min_eff_indices = []
for element, indices in elements.items():
    valid_indices = [i for i in indices if rows[i, 9] not in ["N/A", ""]]
    if valid_indices:
        min_eff_index = min(valid_indices, key=lambda i: float(rows[i, 9]))
        min_eff_indices.append(min_eff_index)

# Find the index of the row with the minimum avg.nu for each element
min_prec_indices = []
for element, indices in elements.items():
    valid_indices = [i for i in indices if rows[i, 10] not in ["N/A", ""]]
    if valid_indices:
        min_prec_index = min(valid_indices, key=lambda i: float(rows[i, 10]))
        min_prec_indices.append(min_prec_index)

for min_eff_index in min_eff_indices:
    for i in range(len(rows[min_eff_index])):
        table[(min_eff_index+1, i)].set_facecolor("lightgreen")

for min_prec_index in min_prec_indices:
    for i in range(len(rows[min_prec_index])):
        table[(min_prec_index+1, i)].set_edgecolor("red")
        table[(min_prec_index+1, i)].set_linewidth(2)

# Add a thick border line as a separator when the element changes
for i in range(1, len(rows)):
    current_element = rows[i, 0].split('.')[0]
    previous_element = rows[i - 1, 0].split('.')[0]
    if current_element != previous_element and previous_element != "":
        for j in range(len(rows[i])):
            # Add a thicker line to the top of the current row to separate elements
            table[(i+1, j)].visible_edges = 'T'  # Only apply the thick line to the top

# Adjust table font size
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.2)

# Create custom legend
legend_patches = [
    Patch(color='orange', label='Precision'),
    Patch(color='lightgreen', label='Efficiency')
]

## Add legend to the plot
#ax.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, 1.15), fontsize=12)

# Display the table
plt.savefig("table.png", bbox_inches='tight', dpi=300)
plt.close()

