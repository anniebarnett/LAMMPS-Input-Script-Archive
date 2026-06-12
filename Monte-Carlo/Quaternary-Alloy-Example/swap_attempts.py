#this script checks how many swap attempts were preformed per species pair. 

import re
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import itertools
import os

# -------------------- Inputs --------------------
log_file = "/home/abarne65/scratch16-mfalk1/annie/defects/4_CrTaVW/run1/log.lammps"

output_dir = os.path.dirname(log_file)

atom_types = [1, 2, 3, 4]  # MoNbTaWV
# atom_types = [1, 3]      # MoTa

species_map = {
    1: "V",
    2: "Cr",
    3: "Ta",
    4: "W"
}

swap_attempts_per_loop = 3750

# -------------------- Parse log --------------------
swap_counts = Counter()
t1, t2 = None, None

with open(log_file, "r") as f:
    for line in f:
        line = line.strip()

        if line.startswith("t1:"):
            t1 = int(line.split(":")[1].strip())

        elif line.startswith("t2:"):
            t2 = int(line.split(":")[1].strip())

            if t1 is not None:
                pair = tuple(sorted((t1, t2)))
                swap_counts[pair] += 1
                t1, t2 = None, None

# -------------------- Prepare pairs --------------------
pairs = sorted(
    tuple(sorted((a, b)))
    for a, b in itertools.combinations_with_replacement(atom_types, 2)
)

attempt_counts = [
    swap_counts.get(pair, 0) * swap_attempts_per_loop
    for pair in pairs
]

labels = [f"{species_map[p[0]]}-{species_map[p[1]]}" for p in pairs]

# Scale factor for plotting
scale = 1e5
scaled_counts = [val / scale for val in attempt_counts]

# -------------------- Plot --------------------
fig, ax = plt.subplots(figsize=(7, 6))

ax.bar(labels, scaled_counts, color="skyblue")

ax.set_xlabel("Atom pair", fontsize=13)
ax.set_ylabel("Total MC Swap Attempts (×10⁵)", fontsize=13)
ax.set_title("MC Swap Attempts per Atom Pair", fontsize=14)
ax.tick_params(axis='x', rotation=45)

# annotate bars
max_val = max(scaled_counts) if scaled_counts else 1

for i, val in enumerate(scaled_counts):
    ax.text(
        i,
        val + max_val * 0.01,
        f"{val:.2f}",
        ha='center',
        va='bottom',
        fontsize=9
    )

plt.tight_layout()

# -------------------- Save + close --------------------
out_path = os.path.join(output_dir, "mc_swap_attempts_per_pair.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved figure:\n{out_path}")
