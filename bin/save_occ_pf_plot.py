#!/usr/bin/env python3
"""Save Occupied% x PassFilter% scatter plot from an Illumina run's InterOp metrics.

Extracted verbatim from ucladx/bcl-qc bclqc.py parse_run_metrics + save_occ_pf_plot
so Nextflow can invoke it directly instead of importing a driver module.
"""
import argparse
import logging
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from interop import py_interop_run_metrics, py_interop_run, py_interop_table
from numpy import zeros, float32

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def parse_run_metrics(run_dir):
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    valid_to_load[py_interop_run.ExtendedTile] = 1
    valid_to_load[py_interop_run.Tile] = 1
    valid_to_load[py_interop_run.Extraction] = 1

    try:
        run_metrics.read(run_dir, valid_to_load)
    except Exception as e:
        logging.error(f"Error occurred trying to open {run_dir}: {e}")
        return None

    columns = py_interop_table.imaging_column_vector()
    py_interop_table.create_imaging_table_columns(run_metrics, columns)
    if columns.size() == 0:
        logging.warning("No interop data found in run metrics.")
        return None

    headers = []
    for i in range(columns.size()):
        column = columns[i]
        if column.has_children():
            headers.extend([f"{column.name()} ({s})" for s in column.subcolumns()])
        else:
            headers.append(column.name())

    column_count = py_interop_table.count_table_columns(columns)
    row_offsets = py_interop_table.map_id_offset()
    py_interop_table.count_table_rows(run_metrics, row_offsets)
    data = zeros((row_offsets.size(), column_count), dtype=float32)
    py_interop_table.populate_imaging_table_data(run_metrics, columns, row_offsets, data.ravel())
    return pd.DataFrame(data, columns=headers)


def save_occ_pf_plot(run_dir, output_dir):
    df = parse_run_metrics(run_dir)
    if df is None:
        logging.warning("Unable to parse Interop files, skipping plot generation.")
        return 1

    os.makedirs(output_dir, exist_ok=True)
    for view in ["Lane"]:
        sns.scatterplot(data=df, x="% Occupied", y="% Pass Filter", hue=view, alpha=0.5, s=8)
        plt.xlim([0, 100])
        plt.ylim([50, 100])
        plt.legend(title=view, bbox_to_anchor=[1.2, 0.9])
        plt.tight_layout()
        out = os.path.join(output_dir, f"occ_pf_{view.lower()}_mqc.jpg")
        plt.savefig(out, dpi=300)
        plt.close()
        logging.info(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    sys.exit(save_occ_pf_plot(**{k.replace("-", "_"): v for k, v in vars(ap.parse_args()).items()}))
