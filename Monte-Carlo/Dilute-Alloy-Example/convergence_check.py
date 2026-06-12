import matplotlib.pyplot as plt
import os

# -------------------- Inputs --------------------
base_dir = "./annie/defects/4_MoNbTaW/"
log_files = [
    os.path.join(base_dir, "run1", "log.lammps"),
    os.path.join(base_dir, "run2", "log.lammps"),
    os.path.join(base_dir, "run3", "log.lammps"),
]

output_dir = base_dir  # change if you want separate folder

# -------------------- Storage --------------------
mc_steps = []
accept_prob = []

energy_steps = []
c_eatoms = []

# -------------------- Parse function --------------------
def parse_log(log_file):
    mc_steps_local = []
    accept_prob_local = []
    energy_steps_local = []
    c_eatoms_local = []

    with open(log_file, "r") as f:
        in_mc_block = False
        headers = []
        idx_step_mc = idx_f1 = idx_f2 = None
        idx_step_en = idx_ce = None

        for line in f:
            line = line.strip()

            # --- MC block detection ---
            if line.startswith("index:"):
                in_mc_block = False

            elif "Step" in line and "f_swaprand[1]" in line:
                headers = line.split()
                idx_step_mc = headers.index("Step")
                idx_f1 = headers.index("f_swaprand[1]")
                idx_f2 = headers.index("f_swaprand[2]")
                in_mc_block = True

            elif in_mc_block:
                if not line or not line[0].isdigit():
                    in_mc_block = False
                    continue
                cols = line.split()
                attempted = float(cols[idx_f1])
                accepted = float(cols[idx_f2])
                if attempted > 0:
                    mc_steps_local.append(1)
                    accept_prob_local.append(accepted / attempted)

            # --- Energy block detection ---
            if "Step" in line and "c_eatoms" in line:
                headers_en = line.split()
                idx_step_en = headers_en.index("Step")
                idx_ce = headers_en.index("c_eatoms")

            elif idx_step_en is not None and idx_ce is not None:
                if not line or not line[0].isdigit():
                    continue
                cols = line.split()
                try:
                    energy_steps_local.append(float(cols[idx_step_en]))
                    c_eatoms_local.append(float(cols[idx_ce]))
                except ValueError:
                    continue

    return mc_steps_local, accept_prob_local, energy_steps_local, c_eatoms_local


# -------------------- Read logs --------------------
mc_counter = 0

for log_file in log_files:
    mc_s, acc_p, en_s, ce = parse_log(log_file)

    mc_steps.extend([mc_counter + i + 1 for i in range(len(mc_s))])
    mc_counter += len(mc_s)
    accept_prob.extend(acc_p)

    energy_steps.extend(en_s)
    c_eatoms.extend(ce)


title_text = r"MoNbTaW, 300°C"

# ==================================================
# 1. ACCEPTANCE PROBABILITY PLOT
# ==================================================
fig1 = plt.figure(figsize=(6, 4))
ax1 = fig1.add_subplot(111)

ax1.plot(mc_steps, accept_prob, color="k", lw=1.5)
ax1.set_ylabel("MC Swap Acceptance Probability", fontsize=13)
ax1.set_xlabel("MC Step Index", fontsize=13)
ax1.set_ylim(0, 1.05)
ax1.grid(True)

ax1.text(0.98, 0.98, title_text,
         transform=ax1.transAxes,
         fontsize=12, ha='right', va='top')

fig1.tight_layout()

out1 = os.path.join(output_dir, "mc_acceptance_probability.png")
fig1.savefig(out1, dpi=300, bbox_inches="tight")
plt.close(fig1)

# ==================================================
# 2. ENERGY CONVERGENCE PLOT
# ==================================================
fig2 = plt.figure(figsize=(6, 4))
ax2 = fig2.add_subplot(111)

ax2.plot(energy_steps, c_eatoms, color="navy", lw=1.5)
ax2.set_xlabel("MD Step", fontsize=13)
ax2.set_ylabel("Potential Energy (eV)", fontsize=13)
ax2.grid(False)

ax2.text(0.98, 0.98, title_text,
         transform=ax2.transAxes,
         fontsize=12, ha='right', va='top')

fig2.tight_layout()

out2 = os.path.join(output_dir, "energy_convergence.png")
fig2.savefig(out2, dpi=300, bbox_inches="tight")
plt.close(fig2)
