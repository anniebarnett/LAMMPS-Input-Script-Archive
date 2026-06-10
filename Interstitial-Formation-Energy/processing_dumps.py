import os
import glob
import re

initial_dir = "initial_dumps"
final_dir = "final_dumps"
output_root = "sorted_runs"

os.makedirs(output_root, exist_ok=True)

pattern = re.compile(r"(.*)_SIA_(\d+)\.dump")

def index_files(directory):
    files = glob.glob(os.path.join(directory, "*.dump"))
    indexed = {}
    for f in files:
        name = os.path.basename(f)
        match = pattern.match(name)
        if match:
            base = match.group(1) + "_SIA_" + match.group(2)
            indexed[base] = os.path.abspath(f)  # absolute path is safer for symlinks
    return indexed

initial_files = index_files(initial_dir)
final_files = index_files(final_dir)

all_keys = sorted(set(initial_files.keys()) & set(final_files.keys()))

for i, key in enumerate(all_keys, start=1):
    run_dir = os.path.join(output_root, f"run{i}")
    os.makedirs(run_dir, exist_ok=True)

    init_src = initial_files[key]
    final_src = final_files[key]

    init_dst = os.path.join(run_dir, "0.dump")
    final_dst = os.path.join(run_dir, "1.dump")

    # remove existing links/files if rerunning
    for dst in [init_dst, final_dst]:
        if os.path.exists(dst):
            os.remove(dst)

    # create symbolic links instead of copying
    os.symlink(init_src, init_dst)
    os.symlink(final_src, final_dst)

print(f"Linked {len(all_keys)} runs into {output_root}/ (no data duplicated)")

