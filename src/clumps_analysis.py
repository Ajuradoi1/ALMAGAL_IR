# Loading libraries

import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable
from matplotlib import pyplot as plt
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# vstack is used to concatenate tables
from astropy.table import vstack
from scipy.stats import gaussian_kde
import sys
import matplotlib

matplotlib.use("Agg")


def source_match(region, alma_sources, spitzer_sources):

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
        ],
        units=[None, None, None, "deg", "deg", "deg", "deg", "arcsec deg"],
        dtype=[
            "str",
            "int64",
            "str",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
        ],
    )

    spitzer_sources_cut = spitzer_sources[
        spitzer_sources["internal_id_01"] == str(region["INTERNAL_ID"])
    ]
    alma_sources_cut = alma_sources[
        alma_sources["CLUMP"] == str(region["INTERNAL_ID"])
    ].reset_index(drop=True)

    match_table_local = match_2_cat(
        cat_1=alma_sources_cut,
        cat_2=spitzer_sources_cut,
        match_table_local=match_table_local,
    )
    return match_table_local


def match_2_cat(cat_1, cat_2, match_table_local):

    if isinstance(cat_1, pd.DataFrame):
        cat_1_id = "CLUMP"
        cat_1_ra = "RA"
        cat_1_dec = "DEC"
    else:
        cat_1_id = "designation"
        cat_1_ra = "ra"
        cat_1_dec = "dec"

    for cat_1_index in range(len(cat_1)):
        cat_1_coords = SkyCoord(
            cat_1[cat_1_ra][cat_1_index],
            cat_1[cat_1_dec][cat_1_index],
            unit=(u.deg, u.deg),
            frame="icrs",
        )
        radius = (
            0.6 * u.arcsec + 0.25 * u.arcsec
        )  # Spitzer resolution + ALMA resolution (assuming 0.25" for ALMA)
        cat_2_coords = SkyCoord(
            cat_2["ra"],
            cat_2["dec"],
            unit=(u.deg, u.deg),
            frame="icrs",
        )
        seps = cat_1_coords.separation(cat_2_coords)

        for cat_2_index in range(len(cat_2)):
            # Only consider the spitzer source with the minimum separation to the alma source

            if seps[cat_2_index] < radius and seps[cat_2_index] == np.min(seps):
                match_table_local.add_row(
                    [
                        cat_1[cat_1_id][cat_1_index],
                        (cat_1_index + 1),
                        cat_2["designation"][cat_2_index],
                        cat_1_coords.ra.deg * u.deg,
                        cat_1_coords.dec.deg * u.deg,
                        cat_2["ra"][cat_2_index] * u.deg,
                        cat_2["dec"][cat_2_index] * u.deg,
                        (seps[cat_2_index] * 3600) * u.arcsec,
                    ]
                )
    return match_table_local


def match_table_fun(alma_regions_file: str, alma_sources_file: str, catalog: str):

    alma_regions = Table.read("../alma_input/" + alma_regions_file, format="ascii.ipac")

    spitzer_sources = Table.read(
        "../output/tables/spitzer_sources_" + catalog + ".tbl", format="ascii.ipac"
    )

    alma_sources = pd.read_csv(
        "../alma_input/" + alma_sources_file,
        comment="\\",
        skiprows=[2],
        delimiter="\s+",
    )

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
        ],
        units=[None, None, None, "deg", "deg", "deg", "deg", "arcsec deg"],
        dtype=[
            "str",
            "int64",
            "str",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
        ],
    )

    test_regions = len(
        alma_regions
    )  # Set this to a smaller number for testing, e.g., 10 or 20
    # test_regions = 100
    regions_subset = alma_regions[:test_regions]

    if os.name == "nt":  # For Windows
        print(
            "Running on Windows. Parallel work is not supported. Processing regions sequentially..."
        )
        for region in tqdm(regions_subset, total=len(regions_subset)):
            result = source_match(region, alma_sources, spitzer_sources)
            match_table = vstack([match_table, result])
    else:
        with Pool(cpu_count()) as pool:
            print(cpu_count(), "cores available. Processing regions in parallel...")

            args = [
                (region, alma_sources, spitzer_sources) for region in regions_subset
            ]

            results = tqdm(
                pool.starmap(source_match, args),
                total=len(regions_subset),
            )
        for result in results:
            match_table = vstack([match_table, result])

    match_table["separation_arcsec"].format = "{:.2E}"

    # print(match_table)

    match_table.write(
        "../output/tables/match_table_" + catalog + ".tbl",
        format="ascii.ipac",
        overwrite=True,
    )
    return match_table


os.makedirs("../output", exist_ok=True)
os.makedirs("../output/figures", exist_ok=True)


def plot_separation_histogram(match_table, bins, catalog):
    plt.figure(figsize=(8, 6))
    plt.hist(
        match_table["separation_arcsec"], bins=bins, color="skyblue", edgecolor="black"
    )
    plt.title("Distribution of Separations Between ALMA and Spitzer Sources")
    plt.xlabel("Separation (arcseconds)")
    plt.ylabel("Number of Matches")
    plt.axvline(
        0.85, color="black", linestyle="--", label="Matching Radius (0.85 arcsec)"
    )
    plt.axvline(0.25, color="red", linestyle="-", label="ALMA Resolution (0.25 arcsec)")
    plt.legend()
    plt.grid(axis="y", alpha=0.75)
    plt.savefig("../output/figures/match_" + catalog + "_hist.svg", format="svg")


def plot_source_matches_all(match_table, catalog):

    _, ax = plt.subplots()

    ax.scatter([0], [0], color="red", label="ALMA Source", alpha=1)

    radius = (
        0.6 * u.arcsec + 0.25 * u.arcsec
    )  # Spitzer resolution + ALMA resolution (assuming 0.25" for ALMA)

    circle_total = plt.Circle(
        (0, 0),
        radius.to_value(u.arcsec),
        color="black",
        fill=False,
        linestyle="--",
        label="Matching Radius (0.85 arcsec)",
    )
    ax.add_patch(circle_total)
    ax.set_ylim(-1, 1)
    ax.set_xlim(1, -1)
    ax.set_ylabel("Declination (arcsec)")
    ax.set_xlabel("Right Ascension (arcsec)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Composition of spitzer sources around ALMA sources")

    for match in match_table:
        ax.scatter(
            (match["spitzer_ra"] - match["alma_ra"])
            * np.cos(np.deg2rad(match["alma_dec"]))
            * 3600,
            (match["spitzer_dec"] - match["alma_dec"]) * 3600,
            color="blue",
            alpha=0.5,
        )
    plt.savefig("../output/figures/match_" + catalog + "_plot.svg", format="svg")


def plot_match_heatmap(match_table, catalog: str):

    def ring(radius, color, linewidth):
        theta = np.linspace(0, 2 * np.pi, 500)
        x_r = radius * np.cos(theta)
        y_r = radius * np.sin(theta)
        plt.plot(x_r, y_r, color=color, linewidth=linewidth, linestyle="--")

    x = []
    y = []
    radius = (
        0.6 * u.arcsec + 0.25 * u.arcsec
    )  # Spitzer resolution + ALMA resolution (assuming 0.25" for ALMA)

    for match in match_table:
        x.append(
            (
                (match["spitzer_ra"] - match["alma_ra"])
                * 3600
                * np.cos(np.deg2rad(match["alma_dec"]))
            )
        )
        y.append(((match["spitzer_dec"] - match["alma_dec"]) * 3600))

    fig, ax = plt.subplots()

    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)

    # grid
    xi, yi = np.mgrid[min(x) : max(x) : 200j, max(y) : min(y) : 200j]
    zi = kde(np.vstack([xi.flatten(), yi.flatten()]))

    im = ax.imshow(
        zi.reshape(xi.shape), origin="lower", extent=[1, -1, -1, 1], cmap="viridis"
    )

    ring(radius, "black", 2)
    ring(0.25, "red", 1)
    ax.scatter(0, 0, color="red", linewidth=1, s=300, marker="+")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(label="Density")

    ax.set_aspect("equal")
    ax.set_ylabel("Declination (arcsec)")
    ax.set_xlabel("Right Ascension (arcsec)")
    ax.set_title("Distribution of Spitzer sources around ALMA sources")

    plt.tight_layout()

    plt.savefig("../output/figures/match_" + catalog + "_heatmap.svg", format="svg")
