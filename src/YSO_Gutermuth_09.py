# YSO classification based on Gutermuth et al. 2009

import pandas as pd
from matplotlib import pyplot as plt
from astropy.table import Table, setdiff
import numpy as np
import os

os.makedirs("../output", exist_ok=True)
os.makedirs("../output/figures", exist_ok=True)


def YSO_Phase_I(spitzer_sources, cat):
    # < 0.2 mag error in all four bands

    good = (  # mask is used to identify null values
        ~spitzer_sources["d3_6m"].mask
        & ~spitzer_sources["d4_5m"].mask
        & ~spitzer_sources["d5_8m"].mask
        & ~spitzer_sources["d8_0m"].mask
    )

    bad = (
        (spitzer_sources["d3_6m"] > 0.2)
        & (spitzer_sources["d4_5m"] > 0.2)
        & (spitzer_sources["d5_8m"] > 0.2)
        & (spitzer_sources["d8_0m"] > 0.2)
    )

    phase_I = spitzer_sources[good & ~bad]

    print(
        f"\n[Phase_I] - {len(phase_I)} sources with good photometry in all four IRAC bands"
    )

    """phase_I.write(
        "../output/tables/phase_I_" + cat + ".tbl",
        format="ascii.ipac",
        overwrite=True,
    )"""

    # Remove PAH galaxies

    _, ax = plt.subplots(2, 1, figsize=(6, 6))

    pah_galaxies = phase_I.copy()

    phase_I = phase_I[
        ~(
            (
                (phase_I["color45_58"] < (1.05 / 1.2) * (phase_I["color58_80"] - 1.0))
                & (phase_I["color45_58"] < 1.05)
                & (phase_I["color58_80"] > 1.0)
                & (phase_I["mag4_5"] > 11.5)
            )
            | (
                (phase_I["color36_58"] < 1.05)
                & (phase_I["mag4_5"] > 11.5)
                & (phase_I["color36_58"] < (1.5 / 2) * (phase_I["color45_80"] - 1.0))
                & (phase_I["color45_80"] > 1.0)
            )
        )
    ]

    pah_galaxies = setdiff(pah_galaxies, phase_I, keys="designation")

    print(f"\n[Phase_I] - {len(pah_galaxies)} PAH galaxies removed")

    _, ax = plt.subplots(2, 1, figsize=(4, 6))

    # Data points

    ax[0].scatter(
        pah_galaxies["color58_80"],
        pah_galaxies["color45_58"],
        s=1,
        color="black",
        alpha=1,
    )

    ax[0].scatter(
        phase_I["color58_80"], phase_I["color45_58"], s=1, color="grey", alpha=0.5
    )

    ax[1].scatter(
        pah_galaxies["color45_80"],
        pah_galaxies["color36_58"],
        s=1,
        color="black",
        alpha=1,
    )

    ax[1].scatter(
        phase_I["color45_80"], phase_I["color36_58"], s=1, color="grey", alpha=0.5
    )

    # Setting up lines to separate PAH galaxies from YSOs

    x = np.linspace(1, 6, 1000)
    cut = 2.2  # This is the x value where the line changes slope
    cut_idx = np.searchsorted(x, cut)  # Find the index of the cut point

    ax[0].plot([1, 1], [-1, 0], color="black", linestyle="-")

    ax[0].plot(
        x[:cut_idx],
        (1.05 / 1.2) * (x[:cut_idx] - 1),
        color="black",
        linestyle="-",
    )
    ax[0].plot(
        x[cut_idx:], np.zeros(len(x[cut_idx:])) + 1.05, color="black", linestyle="-"
    )

    cut = 3  # This is the x value where the line changes slope
    cut_idx = np.searchsorted(x, cut)  # Find the index of the cut point

    ax[1].plot([1, 1], [-2, 0], color="black", linestyle="-")

    ax[1].plot(
        x[:cut_idx],
        (1.5 / 2) * (x[:cut_idx] - 1),
        color="black",
        linestyle="-",
    )
    ax[1].plot(
        x[cut_idx:], np.zeros(len(x[cut_idx:])) + 1.5, color="black", linestyle="-"
    )

    ax[0].set_xlim(-1, 4)
    ax[0].set_ylim(-1, 3)

    ax[1].set_xlim(-2, 6)
    ax[1].set_ylim(-2, 6)

    ax[0].set_xlabel("[5.8] - [8.0]")
    ax[0].set_ylabel("[4.5] - [5.8]")

    ax[1].set_xlabel("[4.5] - [8.0]")
    ax[1].set_ylabel("[3.6] - [5.8]")

    plt.tight_layout()
    plt.savefig(
        "../output/figures/YSO_phase_I_pah_galaxies_" + cat + ".svg", format="svg"
    )

    # Remove AGN sources

    A = 1.95
    A = 0

    agn_sources = phase_I.copy()

    phase_I = phase_I[
        ~(
            (
                (phase_I["color45_80"] > 0.5)
                & (phase_I["mag4_5"] + A > 13.5 + (phase_I["color45_80"] - 2.3) / 0.4)
                & (phase_I["mag4_5"] + A > 13.5)
            )
            & (
                (phase_I["mag4_5"] + A > 14 + (phase_I["color45_80"] - 0.5))
                | (phase_I["mag4_5"] + A > 14.5 - (phase_I["color45_80"] - 1.2) / 0.3)
                | (phase_I["mag4_5"] + A > 14.5)
            )
        )
    ]

    agn_sources = setdiff(agn_sources, phase_I, keys="designation")

    print(f"\n[Phase_I] - {len(agn_sources)} AGN sources removed")

    # Plotting AGN sources

    _, ax = plt.subplots(figsize=(4, 6), layout="tight")

    ax.scatter(
        phase_I["color45_80"], phase_I["mag4_5"] + A, color="grey", alpha=0.5, s=0.5
    )

    ax.scatter(
        agn_sources["color45_80"],
        agn_sources["mag4_5"] + A,
        color="black",
        s=10,
        alpha=1,
    )

    # Plotting line

    x = np.linspace(0.5, 5, 1000)

    cut_1 = 1
    cut_2 = 1.2
    cut_3 = 1.5
    cut_4 = 2.3

    cut_1_idx = np.searchsorted(x, cut_1)
    cut_2_idx = np.searchsorted(x, cut_2)
    cut_3_idx = np.searchsorted(x, cut_3)
    cut_4_idx = np.searchsorted(x, cut_4)

    plt.plot([0.5, 0.5], [18, 14], color="black", linestyle="-")

    plt.plot(x[:cut_1_idx], 14 + x[:cut_1_idx] - 0.5, color="black", linestyle="-")
    plt.plot(
        x[cut_1_idx:cut_2_idx],
        np.zeros(len(x[cut_1_idx:cut_2_idx])) + 14.5,
        color="black",
        linestyle="-",
    )
    plt.plot(
        x[cut_2_idx:cut_3_idx],
        14.5 - (x[cut_2_idx:cut_3_idx] - 1.2) / 0.3,
        color="black",
        linestyle="-",
    )
    plt.plot(
        x[cut_3_idx:cut_4_idx],
        np.zeros(len(x[cut_3_idx:cut_4_idx])) + 13.5,
        color="black",
        linestyle="-",
    )
    plt.plot(
        x[cut_4_idx:], 13.5 + (x[cut_4_idx:] - 2.3) / 0.4, color="black", linestyle="-"
    )

    ax.set_xlabel("[4.5] - [8.0]")
    ax.set_ylabel("[4.5]")
    ax.set_xlim([-2, 6])
    ax.set_ylim([4, 18])
    ax.yaxis.set_inverted(True)

    plt.savefig("../output/figures/YSO_phase_I_AGN_" + cat + ".svg", format="svg")

    # Remove shock emission dominated sources

    phase_I = phase_I[
        (
            ~(
                (
                    phase_I["color36_45"]
                    > (12.0 / 0.55) * (phase_I["color45_80"] - 0.3) + 0.8
                )
                & (phase_I["color45_58"] <= 0.85)
                & (phase_I["color36_45"] > 1.05)
            )
        )
    ]

    # Remove PAH contaminated sources

    sigma_1 = phase_I["d4_5m"] + phase_I["d5_8m"]
    sigma_2 = phase_I["d3_6m"] + phase_I["d4_5m"]

    phase_I = phase_I[
        (
            ~(
                (
                    phase_I["color36_45"] - sigma_2
                    <= 1.4 * (phase_I["color45_58"] + sigma_1 - 0.7)
                )
                & (phase_I["color36_45"] - sigma_2 <= 1.65)
            )
        )
    ]

    # Class I

    phase_I_class_I = phase_I[
        (phase_I["color45_58"] > 0.7) & (phase_I["color36_45"] > 0.7)
    ]
    ids_I = np.array(phase_I_class_I["designation"])

    # Class II
    sigma_2 = phase_I["d3_6m"] + phase_I["d4_5m"]
    sigma_3 = phase_I["d4_5m"] + phase_I["d8_0m"]
    sigma_4 = phase_I["d3_6m"] + phase_I["d5_8m"]

    phase_I_class_II = phase_I[
        (phase_I["color45_80"] - sigma_3 > 0.5)
        & (phase_I["color36_58"] - sigma_4 > 0.35)
        & (
            phase_I["color36_58"] + sigma_4
            <= (0.14 / 0.04) * (phase_I["color45_80"] - sigma_3 - 0.5) + 0.5
        )
        & (phase_I["color36_45"] - sigma_2 > 0.15)
        & ~(np.isin(phase_I["designation"], ids_I))
    ]

    ids_I = np.array(phase_I_class_I["designation"])
    ids_II = np.array(phase_I_class_II["designation"])

    phase_I_no_class = phase_I[
        ~(np.isin(phase_I["designation"], ids_I))
        & ~(np.isin(phase_I["designation"], ids_II))
    ]

    print(
        f"\n[Phase I] - Class I: {len(phase_I_class_I)}, Class II: {len(phase_I_class_II)}, Not Classified: {len(phase_I_no_class)}"
    )

    _, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(
        phase_I["color45_58"], phase_I["color36_45"], s=1, color="black", alpha=0.5
    )
    ax.scatter(
        phase_I_class_I["color45_58"],
        phase_I_class_I["color36_45"],
        s=1,
        color="blue",
    )
    ax.scatter(
        phase_I_class_II["color45_58"],
        phase_I_class_II["color36_45"],
        s=1,
        color="red",
    )

    x_min, x_max = np.min(phase_I["color45_58"]), np.max(phase_I["color45_58"])
    y_min, y_max = np.min(phase_I["color36_45"]), np.max(phase_I["color36_45"])

    x = np.linspace(x_min, x_max, 100)

    ax.axvline(1.05, color="black", linestyle="--")

    ax.set_aspect("equal", adjustable="box")
    ax.set_ylabel("[3.6] - [4.5]")
    ax.set_xlabel("[4.5] - [5.8]")
    ax.set_title("YSO candidates")

    plt.savefig("../output/figures/YSO_phase_I_" + cat + ".png", format="png", dpi=300)

    return phase_I_class_I, phase_I_class_II, phase_I_no_class


def YSO_Phase_II(spitzer_sources, cat):

    # j_msigcom|h_msigcom|k_msigcom|

    phase_II = spitzer_sources[spitzer_sources["has_twomass"] == "True"]

    print(phase_II)

    phase_II = phase_II[
        (phase_II["mag5_8"] == np.nan)
        | (
            phase_II["mag8_0"] == np.nan
        )  # Only consider sources that lack detections either at 5.8 or 8.0 microns
    ]

    phase_II = phase_II[
        (
            (phase_II["h_msigcom"] < 0.1) & (phase_II["k_msigcom"] < 0.1)
        )  # H and K are required at minimum
        & (
            ~(phase_II["j_msigcom"] > 0.1)
        )  # J can be missing, but if present should be good
    ]

    # Preamble done

    # If we have J we use this method

    # line of sight extintion color is

    EJH_HK = 0.11 / 0.065

    HK_color = phase_II["magH"] - phase_II["magK"]
    JH_color = phase_II["magJ"] - phase_II["magH"]

    HK_color_intrinsic = (JH_color - EJH_HK * HK_color - 0.52) / (0.58 - EJH_HK)
    JH_color_intrinsic = np.zeros(len(HK_color_intrinsic))

    for idx, value in enumerate(HK_color_intrinsic):
        if value <= 0.14:
            JH_color_intrinsic[idx] = 0.6
        else:
            JH_color_intrinsic[idx] = 0.58 * value + 0.52

    return 0
