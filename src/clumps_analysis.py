# Loading libraries

import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable, unique
from matplotlib import pyplot as plt
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# vstack is used to concatenate tables
from astropy.table import vstack
from scipy.stats import gaussian_kde
import sys
import matplotlib

matplotlib.use("Agg")

# MARK: ALMA_SOURCE_MATCH

def alma_source_match(region, alma_sources, spitzer_sources):

    match_table_local = QTable(
        names=[
            "region_id",
            "alma_source",
            "spitzer_source",
            "alma_ra",
            "alma_dec",
            "spitzer_ra",
            "spitzer_dec",
            "separation_arcsec",
            "alma_name",
        ],
        units=[None, None, None, "deg", "deg", "deg", "deg", "arcsec deg", None],
        dtype=[
            "str",
            "int64",
            "str",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "str",
        ],
    )

    spitzer_sources_cut = spitzer_sources[
        spitzer_sources["internal_id_01"] == str(region["INTERNAL_ID"])
    ]
    alma_sources_cut = alma_sources[alma_sources["CLUMP"] == str(region["INTERNAL_ID"])]

    alma_coords = SkyCoord(
        alma_sources_cut["RA"],
        alma_sources_cut["DEC"],
        unit=(u.deg, u.deg),
        frame="icrs",
    )
    spitzer_coords = SkyCoord(
        spitzer_sources_cut["ra"],
        spitzer_sources_cut["dec"],
        unit=(u.deg, u.deg),
        frame="icrs",
    )
    radius = 0.6 * u.arcsec  # Spitzer resolution (higher than ALMA's)

    for alma_source, alma_coord in enumerate(alma_coords, start=0):
        seps = alma_coords[alma_source].separation(spitzer_coords)
        for spitzer_source, spitzer_coord in enumerate(spitzer_coords, start=0):
            if seps[spitzer_source] < radius and seps[spitzer_source] == np.min(seps):
                match_table_local.add_row(
                    [
                        alma_sources_cut["CLUMP"][0],
                        alma_sources_cut["ID"][alma_source],
                        spitzer_sources_cut["designation"][spitzer_source],
                        alma_coord.ra.deg * u.deg,
                        alma_coord.dec.deg * u.deg,
                        spitzer_coord.ra.deg * u.deg,
                        spitzer_coord.dec.deg * u.deg,
                        (seps[spitzer_source] * 3600) * u.arcsec,
                        str(alma_sources_cut["CLUMP"][0])
                        + "_"
                        + str(alma_sources_cut["ID"][alma_source]),
                    ]
                )

    return match_table_local

# MARK: ALMA_MATCH_TABLE_FUN

def alma_match_table_fun(spitzer_sources: Table, alma_regions: Table, alma_sources: Table, catalog: str):

    match_table = QTable(
        names=[
            "region_id",
            "alma_source",
            "spitzer_source",
            "alma_ra",
            "alma_dec",
            "spitzer_ra",
            "spitzer_dec",
            "separation_arcsec",
            "alma_name",
        ],
        units=[None, None, None, "deg", "deg", "deg", "deg", "arcsec deg", None],
        dtype=[
            "str",
            "int64",
            "str",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "str",
        ],
    )

    regions_subset = alma_regions

    if os.name == "nt":  # For Windows
        print(
            "Running on Windows. Parallel work is not supported. Processing regions sequentially..."
        )
        for region in tqdm(regions_subset, total=len(regions_subset)):
            result = alma_source_match(region, alma_sources, spitzer_sources)
            match_table = vstack([match_table, result])
    else:
        with Pool(cpu_count()) as pool:
            print(cpu_count(), "cores available. Processing regions in parallel...")

            args = [
                (region, alma_sources, spitzer_sources) for region in regions_subset
            ]

            results = tqdm(
                pool.starmap(alma_source_match, args),
                total=len(regions_subset),
            )
        for result in results:
            match_table = vstack([match_table, result])

    match_table["separation_arcsec"].format = "{:.2E}"

    # print(match_table)

    # Sort by separation first
    match_table.sort("separation_arcsec")

    # Keep first occurrence of each spitzer source
    match_table_unique = unique(match_table, keys="spitzer_source", keep="first")

    match_table_unique.write(
        "../output/tables/match_table_" + catalog + ".tbl",
        format="ascii.ipac",
        overwrite=True,
    )
    return match_table_unique


os.makedirs("../output", exist_ok=True)
os.makedirs("../output/figures", exist_ok=True)

# MARK: ALMA_MATCH_PLOTS

def plot_separation_histogram(match_table, bin_width, catalog):
    plt.figure(figsize=(5, 4))
    plt.hist(
        match_table["separation_arcsec"],
        bins=np.array(range(0, 65, bin_width)) / 100,
        color="skyblue",
        edgecolor="black",
    )
    # plt.title("Distribution of Separations Between ALMA and Spitzer Sources")
    plt.xlabel("Separation (arcseconds)")
    plt.ylabel("Matches")
    plt.axvline(
        0.6, color="black", linestyle="--", label="Matching Radius (0.60 arcsec)"
    )
    plt.axvline(0.25, color="red", linestyle="-", label="ALMA Resolution (0.25 arcsec)")

    plt.grid(axis="y", alpha=0.75)
    plt.savefig("../output/figures/"+catalog+"/match_" + catalog + "_hist.pdf", format="pdf")


def plot_source_matches_all(match_table, catalog):

    _, ax = plt.subplots()

    ax.scatter([0], [0], color="red", label="ALMA Source", alpha=1)

    radius = 0.6 * u.arcsec  # Spitzer resolution

    circle_total = plt.Circle(
        (0, 0),
        radius.to_value(u.arcsec),
        color="black",
        fill=False,
        linestyle="--",
        label="Matching Radius (0.6 arcsec)",
    )
    ax.add_patch(circle_total)
    ax.set_ylim(-1, 1)
    ax.set_xlim(1, -1)
    ax.set_ylabel("Declination (arcsec)")
    ax.set_xlabel("Right Ascension (arcsec)")
    ax.set_aspect("equal", adjustable="box")
    # ax.set_title("Composition of spitzer sources around ALMA sources")

    for match in match_table:
        ax.scatter(
            (match["spitzer_ra"] - match["alma_ra"])
            * np.cos(np.deg2rad(match["alma_dec"]))
            * 3600,
            (match["spitzer_dec"] - match["alma_dec"]) * 3600,
            color="blue",
            alpha=0.5,
        )
    plt.savefig("../output/figures/"+catalog+"/match_" + catalog + "_plot.pdf", format="pdf")


def plot_match_heatmap(match_table, catalog: str):

    def ring(radius, color, linewidth, alpha):
        theta = np.linspace(0, 2 * np.pi, 500)
        x_r = radius * np.cos(theta)
        y_r = radius * np.sin(theta)
        plt.plot(
            x_r, y_r, color=color, linewidth=linewidth, linestyle="--", alpha=alpha
        )

    x = []
    y = []
    radius = 0.6 * u.arcsec  # Spitzer resolution
    fig, ax = plt.subplots()

    for match in match_table:
        x.append(
            (
                (match["spitzer_ra"] - match["alma_ra"])
                * 3600
                * np.cos(np.deg2rad(match["alma_dec"]))
            )
        )
        y.append(((match["spitzer_dec"] - match["alma_dec"]) * 3600))

    # x, y = your data

    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)

    # fixed grid: x in [-1,1], y in [-1,1]
    xi, yi = np.mgrid[-0.65:0.65:200j, -0.65:0.65:200j]

    zi = kde(np.vstack([xi.ravel(), yi.ravel()])).reshape(xi.shape)

    fig, ax = plt.subplots(tight_layout=True)

    im = ax.imshow(
        zi.T,  # transpose is essential
        origin="lower",
        extent=[-0.65, 0.65, -0.65, 0.65],
        cmap="viridis",
        aspect="equal",
    )

    # overlay points
    ax.scatter(x, y, color="black", alpha=0.6, s=1,marker='x',linewidth=10)

    # flip x so it goes 1 -> -1
    ax.set_xlim(0.65, -0.65)
    ax.set_ylim(-0.65, 0.65)

    ring(radius, "white", linewidth=2, alpha=1)
    ring(0.25, "red", linewidth=2, alpha=0.8)
    ax.scatter(0, 0, color="red", linewidth=1, s=300, marker="+")

    ax.set_aspect("equal", adjustable="box")
    ax.set_ylabel("Declination (arcsec)")
    ax.set_xlabel("Right Ascension (arcsec)")
    # ax.set_title("Distribution of Spitzer sources around ALMA sources")

    fig.colorbar(im, ax=ax, label="Density")
    fig.savefig("../output/figures/"+catalog+"/match_" + catalog + "_heatmap.pdf", format="pdf")
