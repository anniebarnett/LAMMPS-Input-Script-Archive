from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
import numpy as np
from ovito.data import DataCollection
import os
from multiprocessing import Pool
from tqdm.notebook import tqdm
import time
import pandas as pd
import glob
from multiprocessing import Pool, Manager
from filelock import FileLock

import pyarrow as pa
import pyarrow.parquet as pq

base_dir = '/home/abarne65/scratch16-mfalk1/annie/defects/dumbbells/dilute/Mo/sorted_runs/'
output_parquet = os.path.join(base_dir, 'D_Mo_DBS.parquet')

run_dirs = sorted(glob.glob(os.path.join(base_dir, 'run*')))

all_rows = []

for input_folder in run_dirs:

    run_name = os.path.basename(input_folder)
    run_number = int(run_name.replace('run',''))

    dump_files = sorted(glob.glob(os.path.join(input_folder, '*.dump')))
    if not dump_files:
        continue

    pipeline = import_file(dump_files, multiple_frames=True)

    pipeline.modifiers.append(
        WignerSeitzAnalysisModifier(
            output_displaced=True,
            reference_frame=0,
            affine_mapping=ReferenceConfigurationModifier.AffineMapping.ToReference
        )
    )

    # ✅ ONLY FINAL FRAME
    frame = pipeline.source.num_frames - 1
    data = pipeline.compute(frame)

    timestep = data.attributes["Timestep"]

    # Pull arrays ONCE
    occ = np.asarray(data.particles['Occupancy'])

    # ✅ mask BEFORE building dataframe
    mask = occ > 1

    if not np.any(mask):
        continue

    pos = np.asarray(data.particles['Position'])[mask]
    ptype = np.asarray(data.particles['Particle Type'])[mask]

    df = pd.DataFrame({
        'run': run_number,
        'timestep': timestep,
        'Particle_Type': ptype,
        'Position_X': pos[:,0],
        'Position_Y': pos[:,1],
        'Position_Z': pos[:,2],
        'Occupancy': occ[mask]
    })

    all_rows.append(df)

    print(f"{run_name}: {len(df)} interstitials")

# Combine and write
if all_rows:
    final_df = pd.concat(all_rows, ignore_index=True)
    final_df = final_df.sort_values(by=['run']).reset_index(drop=True)

    pq.write_table(pa.Table.from_pandas(final_df), output_parquet, compression='snappy')
    print(f"\nSaved: {output_parquet}")
