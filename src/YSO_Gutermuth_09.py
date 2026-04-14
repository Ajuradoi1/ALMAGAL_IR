# YSO classification based on Gutermuth et al. 2009

import pandas as pd
from matplotlib import pyplot as plt
from astropy.table import Table
import numpy as np
import os

os.makedirs("../output", exist_ok=True)
os.makedirs("../output/figures", exist_ok=True)


def YSO_Phase_I(cat):

    spitzer_sources = Table.read(
        "../output/tables/spitzer_sources_" + cat + ".tbl", format="ipac"
    )

    # < 0.2 mag error in all four bands

    phase_I = spitzer_sources[
        (
            (spitzer_sources["d3_6m"] < 0.2)
            & (spitzer_sources["d4_5m"] < 0.2)
            & (spitzer_sources["d5_8m"] < 0.2)
            & (spitzer_sources["d8_0m"] < 0.2)
        )
    ]

    # Remove PAH galaxies

    phase_I = phase_I[(~((phase_I["color45_58"] < 1.05) & (phase_I["mag4_5"] > 11.5)))]

    phase_I = phase_I[
        (
            ~(
                (phase_I["color45_58"] < (1.05 / 1.2) * (phase_I["color58_80"] - 1.0))
                & (phase_I["color58_80"] > 1.0)
            )
        )
    ]

    phase_I = phase_I[
        (
            ~(
                (phase_I["color36_58"] < (1.5 / 2) * (phase_I["color45_80"] - 1.0))
                & (phase_I["color45_80"] > 1.0)
            )
        )
    ]

    # Remove AGN sources

    phase_I = phase_I[
        (
            ~(
                (
                    (phase_I["color45_80"] > 0.5)
                    & (phase_I["mag4_5"] > 13.5 + (phase_I["color45_80"] - 2.3) / 0.4)
                    & (phase_I["mag4_5"] > 13.5)
                )
                & (
                    (phase_I["mag4_5"] > 14 + (phase_I["color45_80"] - 0.5))
                    | (phase_I["mag4_5"] > 14.5 - (phase_I["color45_80"] - 1.2) / 0.3)
                    | (phase_I["mag4_5"] > 14.5)
                )
            )
        )
    ]

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
        phase_I["color45_58"], phase_I["color36_45"], s=0.8, color="black", alpha=0.5
    )
    ax.scatter(
        phase_I_class_I["color45_58"],
        phase_I_class_I["color36_45"],
        s=0.8,
        color="blue",
    )
    ax.scatter(
        phase_I_class_II["color45_58"],
        phase_I_class_II["color36_45"],
        s=0.8,
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


def YSO_Phase_II(cat):
    return 0
