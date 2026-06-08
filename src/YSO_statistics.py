import pandas as pd
from matplotlib import pyplot as plt
from astropy.table import Table, setdiff
import numpy as np
import os


def YSO_plots(spitzer_sources, alma_proceeding, cat):

    # pah_galaxies|agn_sources|sed_sources|pah_contaminated|p1_class_I|p1_class_II|p1_noclass|
    # p2_class_I|p2_class_II|p2_noclass|

    spitzer_sources.write(
        "../output/tables/WORKING.tbl", format="ascii.ipac", overwrite=True
    )

    # MARK: COUNT

    fig, ax = plt.subplots(tight_layout=True)

    cols = [
        "pah_galaxies",
        "agn_sources",
        "sed_sources",
        "pah_contaminated",
        "p1_class_I",
        "p1_class_II",
        "p1_noclass",
    ]

    if alma_proceeding:
        cols.extend(["p2_class_I", "p2_class_II", "p2_noclass"])
        # alma_sources = Table.read(
        #    format="ascii.ipac",
        # )
        match_table = Table.read(
            "../output/tables/match_table_" + cat.split("_ALMA")[0] + ".tbl",
            format="ascii.ipac",
        )

    x = []
    for element in cols:
        x.append(len(spitzer_sources[spitzer_sources[element]]))
        print(f"{element}: {x[-1]}")

    ax.legend()

    fig.savefig("../output/figures/WORKING.svg", format="svg")
