# YSO classification based on Gutermuth et al. 2009

import pandas as pd
from matplotlib import pyplot as plt
from astropy.table import Table, setdiff
import numpy as np
import os

os.makedirs("../output", exist_ok=True)
os.makedirs("../output/figures", exist_ok=True)

# MARK: PHASE 1


def YSO_Phase_I(spitzer_sources, cat):

    def safe_setdiff(t1, t2, keys):
        if len(t1) == 0:
            return t1
        if len(t2) == 0:
            return t1
        return setdiff(t1, t2, keys=keys)

    # < 0.2 mag error in all four bands

    A3_6 = 0.63
    A4_5 = 0.57
    A5_8 = 0.49
    A8_0 = 0.55

    A3_6, A4_5, A5_8, A8_0 = [0, 0, 0, 0]

    # test purposes

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
        f"\n\n[Phase I] - {len(phase_I)} sources with good photometry in all four IRAC bands"
    )

    """phase_I.write(
        "../output/tables/phase_I_" + cat + ".tbl",
        format="ascii.ipac",
        overwrite=True,
    )"""

    # MARK: PAH GALAXIES
    # Remove PAH galaxies

    # _, ax = plt.subplots(2, 1, figsize=(6, 6))

    pah_galaxies = phase_I.copy()

    phase_I = phase_I[
        ~(
            (
                (
                    phase_I["color45_58"] + (A4_5 - A5_8) * phase_I["Ak_NODEC"]
                    < (1.05 / 1.2)
                    * (
                        phase_I["color58_80"]
                        + (A5_8 - A8_0) * phase_I["Ak_NODEC"]
                        - 1.0
                    )
                )
                & (phase_I["color45_58"] + (A4_5 - A5_8) * phase_I["Ak_NODEC"] < 1.05)
                & (phase_I["color58_80"] + (A5_8 - A8_0) * phase_I["Ak_NODEC"] > 1.0)
                & (phase_I["mag4_5"] + A4_5 * phase_I["Ak_NODEC"] > 11.5)
            )
            | (
                (phase_I["color36_58"] + (A3_6 - A5_8) * phase_I["Ak_NODEC"] < 1.05)
                & (phase_I["mag4_5"] + A4_5 * phase_I["Ak_NODEC"] > 11.5)
                & (
                    phase_I["color36_58"] + (A3_6 - A5_8) * phase_I["Ak_NODEC"]
                    < (1.5 / 2)
                    * (
                        phase_I["color45_80"]
                        + (A4_5 - A8_0) * phase_I["Ak_NODEC"]
                        - 1.0
                    )
                )
                & (phase_I["color45_80"] + (A4_5 - A8_0) * phase_I["Ak_NODEC"] > 1.0)
            )
        )
    ]

    pah_galaxies = safe_setdiff(pah_galaxies, phase_I, keys="designation")

    print(f"\n[Phase I] - {len(pah_galaxies)} PAH galaxies removed")

    spitzer_sources["pah_galaxies"] = False
    spitzer_sources["pah_galaxies"][
        np.isin(spitzer_sources["designation"], pah_galaxies["designation"])
    ] = True

    _, ax = plt.subplots(2, 1, figsize=(4, 6))

    # Data points

    ax[0].scatter(
        pah_galaxies["color58_80"] + (A5_8 - A8_0) * pah_galaxies["Ak_NODEC"],
        pah_galaxies["color45_58"] + (A4_5 - A5_8) * pah_galaxies["Ak_NODEC"],
        s=3,
        color="black",
        alpha=1,
    )

    ax[0].scatter(
        phase_I["color58_80"] + (A5_8 - A8_0) * phase_I["Ak_NODEC"],
        phase_I["color45_58"] + (A4_5 - A5_8) * phase_I["Ak_NODEC"],
        s=3,
        color="grey",
        alpha=0.5,
    )

    ax[1].scatter(
        pah_galaxies["color45_80"] + (A4_5 - A8_0) * pah_galaxies["Ak_NODEC"],
        pah_galaxies["color36_58"] + (A3_6 - A5_8) * pah_galaxies["Ak_NODEC"],
        s=3,
        color="black",
        alpha=1,
    )

    ax[1].scatter(
        phase_I["color45_80"] + (A4_5 - A8_0) * phase_I["Ak_NODEC"],
        phase_I["color36_58"] + (A3_6 - A5_8) * phase_I["Ak_NODEC"],
        s=3,
        color="grey",
        alpha=0.5,
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

    # MARK: AGN

    # Remove AGN sources

    agn_sources = phase_I.copy()

    phase_I = phase_I[
        ~(
            (
                (phase_I["color45_80"] + (A4_5 - A8_0) * phase_I["Ak_NODEC"] > 0.5)
                & (
                    phase_I["mag4_5"] + (A4_5) * phase_I["Ak_NODEC"]
                    > 13.5
                    + (
                        phase_I["color45_80"]
                        + (A4_5 - A8_0) * phase_I["Ak_NODEC"]
                        - 2.3
                    )
                    / 0.4
                )
                & (phase_I["mag4_5"] + (A4_5) * phase_I["Ak_NODEC"] > 13.5)
            )
            & (
                (
                    phase_I["mag4_5"] + A4_5 * phase_I["Ak_NODEC"]
                    > 14
                    + (
                        phase_I["color45_80"]
                        + (A4_5 - A8_0) * phase_I["Ak_NODEC"]
                        - 0.5
                    )
                )
                | (
                    phase_I["mag4_5"] + A4_5 * phase_I["Ak_NODEC"]
                    > 14.5
                    - (
                        phase_I["color45_80"]
                        + (A4_5 - A8_0) * phase_I["Ak_NODEC"]
                        - 1.2
                    )
                    / 0.3
                )
                | (phase_I["mag4_5"] + A4_5 * phase_I["Ak_NODEC"] > 14.5)
            )
        )
    ]

    agn_sources = safe_setdiff(agn_sources, phase_I, keys="designation")

    print(f"\n[Phase I] - {len(agn_sources)} AGN sources removed")

    spitzer_sources["agn_sources"] = False

    spitzer_sources["agn_sources"][
        np.isin(spitzer_sources["designation"], agn_sources["designation"])
    ] = True

    # Plotting AGN sources

    _, ax = plt.subplots(figsize=(4, 6), layout="tight")

    ax.scatter(
        phase_I["color45_80"] + (A4_5 - A8_0) * phase_I["Ak_NODEC"],
        phase_I["mag4_5"] + A4_5 * phase_I["Ak_NODEC"],
        color="grey",
        alpha=0.5,
        s=2,
    )

    ax.scatter(
        agn_sources["color45_80"] + (A4_5 - A8_0) * agn_sources["Ak_NODEC"],
        agn_sources["mag4_5"] + A4_5 * agn_sources["Ak_NODEC"],
        color="black",
        s=3,
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
    plt.tight_layout()
    plt.savefig("../output/figures/YSO_phase_I_AGN_" + cat + ".svg", format="svg")

    # MARK: SED SOURCES

    # Remove shock emission dominated sources

    sed = phase_I.copy()

    phase_I = phase_I[
        (
            ~(
                (phase_I["color45_58"] + (A4_5 - A5_8) * phase_I["Ak_NODEC"] <= 0.85)
                & (phase_I["color36_45"] + (A3_6 - A4_5) * phase_I["Ak_NODEC"] > 1.05)
                & (
                    phase_I["color36_45"] + (A3_6 - A4_5) * phase_I["Ak_NODEC"]
                    > (1.2 / 0.55)
                    * (
                        phase_I["color45_58"]
                        + (A4_5 - A5_8) * phase_I["Ak_NODEC"]
                        - 0.3
                    )
                    + 0.8
                )
            )
        )
    ]

    sed = safe_setdiff(sed, phase_I, keys="designation")
    print(f"\n[Phase I] - {len(sed)} Shock Emission Dominated sources removed")

    spitzer_sources["sed_sources"] = False
    spitzer_sources["sed_sources"][
        np.isin(spitzer_sources["designation"], sed["designation"])
    ] = True

    # MARK: PAH CONTAMINATED

    # Remove PAH contaminated sources

    sigma_1 = phase_I["d4_5m"] + phase_I["d5_8m"]
    sigma_2 = phase_I["d3_6m"] + phase_I["d4_5m"]

    pah_contam = phase_I.copy()

    phase_I = phase_I[
        (
            ~(
                (
                    phase_I["color36_45"]
                    + (A3_6 - A4_5) * phase_I["Ak_NODEC"]
                    - sigma_2
                    <= 1.4
                    * (
                        (phase_I["color45_58"] + (A4_5 - A5_8) * phase_I["Ak_NODEC"])
                        + sigma_1
                        - 0.7
                    )
                )
                & (
                    phase_I["color36_45"]
                    + (A3_6 - A4_5) * phase_I["Ak_NODEC"]
                    - sigma_2
                    <= 1.65
                )
            )
        )
    ]

    pah_contam = safe_setdiff(pah_contam, phase_I, keys="designation")
    print(f"\n[Phase I] - {len(pah_contam)} PAH-contaminated apertures removed")

    spitzer_sources["pah_contaminated"] = False
    spitzer_sources["pah_contaminated"][
        np.isin(spitzer_sources["designation"], pah_contam["designation"])
    ] = True

    # MARK: CLASS I
    # Class I

    class_I = phase_I.copy()

    phase_I = phase_I[
        ~(
            (phase_I["color45_58"] + (A4_5 - A5_8) * phase_I["Ak_NODEC"] > 0.7)
            & (phase_I["color36_45"] + (A3_6 - A4_5) * phase_I["Ak_NODEC"] > 0.7)
        )
    ]

    class_I = safe_setdiff(class_I, phase_I, keys="designation")
    print(f"\n[Phase I] - {len(class_I)} Class I sources")

    # Plotting

    _, ax = plt.subplots(2, 1, figsize=(4, 6), constrained_layout=True)

    ax[0].scatter(
        phase_I["color45_58"] + (A4_5 - A5_8) * phase_I["Ak_NODEC"],
        phase_I["color36_45"] + (A3_6 - A4_5) * phase_I["Ak_NODEC"],
        color="grey",
        alpha=0.5,
        s=2,
    )

    ax[0].scatter(
        class_I["color45_58"] + (A4_5 - A5_8) * class_I["Ak_NODEC"],
        class_I["color36_45"] + (A3_6 - A4_5) * class_I["Ak_NODEC"],
        color="grey",
        alpha=0.5,
        s=2,
    )

    ax[0].scatter(
        sed["color45_58"] + (A4_5 - A5_8) * sed["Ak_NODEC"],
        sed["color36_45"] + (A3_6 - A4_5) * sed["Ak_NODEC"],
        color="black",
        s=3,
        alpha=1,
    )

    ax[0].scatter(
        pah_contam["color45_58"] + (A4_5 - A5_8) * pah_contam["Ak_NODEC"],
        pah_contam["color36_45"] + (A3_6 - A4_5) * pah_contam["Ak_NODEC"],
        color="black",
        s=3,
        alpha=1,
    )

    x = np.linspace(0.42, 0.85, 100)

    # sed line

    ax[0].plot([-1, 0.42], [1.05, 1.05], color="black", linestyle="-")
    ax[0].plot(x, (1.2 / 0.55) * (x - 0.3) + 0.8, color="black", linestyle="-")
    ax[0].plot([0.85, 0.85], [2, 3.5], color="black", linestyle="-")

    # PAH-contaminated line

    x = np.linspace(0.1, 1.77, 100)

    ax[0].plot(x, 1.4 * (x - 0.7) + 0.15, color="black", linestyle="-")
    ax[0].plot([1.77, 3], [1.65, 1.65], color="black", linestyle="-")

    ax[0].set_xlim([-1, 3])
    ax[0].set_ylim([-0.5, 3.5])

    ax[0].set_xlabel("[4.5] - [5.8]")
    ax[0].set_ylabel("[3.6] - [4.5]")

    ax[0].text(-0.5, 1.5, "Shock\nEmission")
    ax[0].text(1, -0.3, "PAH Aper. Contam.")
    ax[1].text(2.3, 2.5, "Class I")

    # Fig b

    ax[1].scatter(
        phase_I["color45_58"] + (A4_5 - A5_8) * phase_I["Ak_NODEC"],
        phase_I["color36_45"] + (A3_6 - A4_5) * phase_I["Ak_NODEC"],
        color="grey",
        alpha=0.5,
        s=2,
    )

    ax[1].scatter(
        class_I["color45_58"] + (A4_5 - A5_8) * class_I["Ak_NODEC"],
        class_I["color36_45"] + (A3_6 - A4_5) * class_I["Ak_NODEC"],
        color="black",
        alpha=1,
        s=3,
    )

    ax[1].scatter(
        sed["color45_58"] + (A4_5 - A5_8) * sed["Ak_NODEC"],
        sed["color36_45"] + (A3_6 - A4_5) * sed["Ak_NODEC"],
        color="grey",
        s=0.2,
        alpha=0.5,
    )

    ax[1].scatter(
        pah_contam["color45_58"] + (A4_5 - A5_8) * pah_contam["Ak_NODEC"],
        pah_contam["color36_45"] + (A3_6 - A4_5) * pah_contam["Ak_NODEC"],
        color="grey",
        s=0.2,
        alpha=0.5,
    )

    # sed line

    x = np.linspace(0.42, 0.85, 100)
    ax[1].plot([-1, 0.42], [1.05, 1.05], color="black", linestyle="-")
    ax[1].plot(x, (1.2 / 0.55) * (x - 0.3) + 0.8, color="black", linestyle="-")
    ax[1].plot([0.85, 0.85], [2, 3.5], color="black", linestyle="-")

    # PAH-contaminated line

    x = np.linspace(0.1, 1.77, 100)
    ax[1].plot(x, 1.4 * (x - 0.7) + 0.15, color="black", linestyle="-")
    ax[1].plot([1.77, 3], [1.65, 1.65], color="black", linestyle="-")

    # Class I selection line

    ax[1].plot([0.7, 1.09], [0.7, 0.7], color="black", linestyle="-")
    ax[1].plot([0.7, 0.7], [0.7, 1.67], color="black", linestyle="-")

    ax[1].set_xlim([-1, 3])
    ax[1].set_ylim([-0.5, 3.5])

    ax[1].set_xlabel("[4.5] - [5.8]")
    ax[1].set_ylabel("[3.6] - [4.5]")

    plt.savefig("../output/figures/YSO_phase_I_SED_PAH_" + cat + ".svg", format="svg")

    # MARK: Class II

    # Class II
    sigma_2 = phase_I["d3_6m"] + phase_I["d4_5m"]
    sigma_3 = phase_I["d4_5m"] + phase_I["d8_0m"]
    sigma_4 = phase_I["d3_6m"] + phase_I["d5_8m"]

    class_II = phase_I.copy()

    phase_I = phase_I[
        ~(
            (
                phase_I["color45_80"] + (A4_5 - A8_0) * phase_I["Ak_NODEC"] - sigma_3
                > 0.5
            )
            & (
                phase_I["color36_58"] + (A3_6 - A5_8) * phase_I["Ak_NODEC"] - sigma_4
                > 0.35
            )
            & (
                phase_I["color36_58"] + (A3_6 - A5_8) * phase_I["Ak_NODEC"] + sigma_4
                <= (0.14 / 0.04)
                * (
                    phase_I["color45_80"]
                    + (A4_5 - A8_0) * phase_I["Ak_NODEC"]
                    - sigma_3
                    - 0.5
                )
                + 0.5
            )
            & (
                phase_I["color36_45"] + (A3_6 - A4_5) * phase_I["Ak_NODEC"] - sigma_2
                > 0.15
            )
        )
    ]

    class_II = safe_setdiff(class_II, phase_I, keys="designation")

    print(f"\n[Phase I] - {len(class_II)} Class II sources")
    print(
        f"\n[Phase I] - Class I: {len(class_I)}, Class II: {len(class_II)}, Not Classified: {len(phase_I)}"
    )

    # Plotting Class II

    _, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)

    ax.scatter(
        phase_I["color45_80"] + (A4_5 - A8_0) * phase_I["Ak_NODEC"],
        phase_I["color36_58"] + (A3_6 - A5_8) * phase_I["Ak_NODEC"],
        s=2,
        alpha=0.5,
        color="grey",
    )
    ax.scatter(
        class_I["color45_80"] + (A4_5 - A8_0) * class_I["Ak_NODEC"],
        class_I["color36_58"] + (A3_6 - A5_8) * class_I["Ak_NODEC"],
        s=0.1,
        alpha=0.5,
        marker="+",
        color="black",
        label="class I",
    )
    ax.scatter(
        class_II["color45_80"] + (A4_5 - A8_0) * class_II["Ak_NODEC"],
        class_II["color36_58"] + (A3_6 - A5_8) * class_II["Ak_NODEC"],
        s=3,
        alpha=1,
        color="black",
        label="class II",
    )

    x = np.linspace(0.5, 6, 100)

    ax.plot([0.5, 6], [0.5, 0.5], color="black", linestyle="-")
    ax.plot(x, (0.14 / 0.04) * (x - 0.5) + 0.5, color="black", linestyle="-")

    ax.set_ylabel("[3.6] - [5.8]")
    ax.set_xlabel("[4.5] - [8.0]")

    ax.set_xlim([-2, 6])
    ax.set_ylim([-2, 6])

    ax.text(4, 1, "Class II")

    plt.savefig("../output/figures/YSO_phase_I_class_II_" + cat + ".svg", format="svg")

    # return class_I, class_II, phase_I

    spitzer_sources["p1_class_I"] = False
    spitzer_sources["p1_class_I"][
        np.isin(spitzer_sources["designation"], class_I["designation"])
    ] = True

    spitzer_sources["p1_class_II"] = False

    spitzer_sources["p1_class_II"][
        np.isin(spitzer_sources["designation"], class_II["designation"])
    ] = True

    spitzer_sources["p1_noclass"] = False

    spitzer_sources["p1_noclass"][
        np.isin(spitzer_sources["designation"], phase_I["designation"])
    ] = True

    return spitzer_sources


# MARK: PHASE II


def YSO_Phase_II(spitzer_sources, cat):

    def safe_setdiff(t1, t2, keys):
        if len(t1) == 0:
            return t1
        if len(t2) == 0:
            return t1
        return setdiff(t1, t2, keys=keys)

    phase_II = spitzer_sources.copy()

    good = (
        (~phase_II["mag5_8"].mask | ~phase_II["mag8_0"].mask)
        & ~np.isnan(phase_II["h_msigcom"])
        & ~np.isnan(phase_II["k_msigcom"])
        & ~phase_II["color36_45"].mask
    )

    bad = (phase_II["h_msigcom"] > 0.1) & (phase_II["k_msigcom"] > 0.1)

    notJ = ~(np.isnan(phase_II["j_msigcom"]))
    phase_II = phase_II[good & ~bad]

    no_class = phase_II.copy()

    no_class.write(
        "../output/tables/phase_II_" + cat + ".tbl", format="ascii.ipac", overwrite=True
    )

    EJH_HK = 1.73

    EK36_HK = 1.0 / 1.49

    EK45_HK = 1.0 / 1.17

    E3645_HK = EK45_HK - EK36_HK

    print(f"[Phase II] - {len(phase_II)} objects to evaluate")

    data_phase_II = {}

    # FAILED CODE to be fixed after Guthermuth 09 approximations are understood

    """for idx, source in enumerate(phase_II):
        j_sig = source["j_msigcom"]
        jh_meas = source["j_m"] - source["h_m"]
        hk_meas = source["h_m"] - source["k_m"]

        if (not np.isnan(j_sig)) and (j_sig < 0.1):
            colorHK_meas = hk_meas
            colorJH_meas = jh_meas

            colorHK_0 = (colorJH_meas - EJH_HK * colorHK_meas - 0.52) / (0.58 - EJH_HK)

            if colorHK_0 <= 0.14:
                colorJH_0 = 0.6
            else:
                colorJH_0 = 0.58 * colorHK_0 + 0.52

        else:
            # still define measured quantities to avoid undefined vars
            colorHK_meas = hk_meas

            colorHK_0 = (
                1.33 * (E3645_HK * colorHK_meas - source["color36_45"]) - 0.133
            ) / (1.33 * E3645_HK - 1)

            colorJH_0 = np.nan  # better than -999 for scientific pipelines

        colorK_36_0 = (
            source["k_m"] - source["mag3_6"] - (colorHK_meas - colorHK_0) * EK36_HK
        )

        color36_45_0 = source["color36_45"] - (colorHK_meas - colorHK_0) * E3645_HK

        if colorJH_0 == np.nan:
            if color36_45_0 <= 0.06:
                colorHK_0 = 0.2
            else:
                colorHK_0 = 1.33 * color36_45_0 + 0.133

        data_phase_II[idx] = [colorHK_0, colorJH_0, colorK_36_0, color36_45_0]
    """

    # WORKAROUND, ONLY ALMA DATA

    # DE-REDDENING OF COLORS [H - K], [K - 3.6], [3.6 - 4.5],

    # Ak_NODEC

    # ratio respect Ak https://arxiv.org/pdf/0910.3037

    A_H = 1.56
    A_K = 1
    A_36 = 0.63
    A_45 = 0.57

    for idx, source in enumerate(phase_II):
        colorHK_0 = source["h_m"] - source["k_m"] + (A_H - A_K) * source["Ak_NODEC"]
        colorK_36_0 = (
            source["k_m"] - source["mag3_6"] + (A_K - A_36) * source["Ak_NODEC"]
        )
        color36_45_0 = source["color36_45"] + (A_36 - A_45) * source["Ak_NODEC"]

        data_phase_II[idx] = [colorHK_0, colorK_36_0, color36_45_0]

    for idx, col in enumerate(
        #    ["colorHK_0", "colorJH_0", "colorK_36_0", "color36_45_0"]
        ["colorHK_0", "colorK_36_0", "color36_45_0"]
    ):
        phase_II[col] = [data_phase_II[i][idx] for i in range(len(phase_II))]

    phase_II.write(
        "../output/tables/TestPhase_II.tbl", format="ascii.ipac", overwrite=True
    )

    _, ax = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

    ax[0].scatter(
        phase_II["h_m"] - phase_II["k_m"],
        phase_II["j_m"] - phase_II["h_m"],
        s=2,
        color="grey",
    )
    ax[0].set_ylabel("J-H")
    ax[0].set_xlabel("H-K")
    ax[0].set_xlim([-0.6, 2.5])
    ax[0].set_ylim([0, 3])

    plt.show()

    # We now apply conditions

    sigma_1 = phase_II["d3_6m"] + phase_II["d4_5m"]

    sigma_2 = phase_II["k_msigcom"] + phase_II["d3_6m"]

    class_II = phase_II.copy()

    phase_II = phase_II[
        ~(
            (phase_II["color36_45_0"] - sigma_1 > 0.101)
            & (phase_II["colorK_36_0"] - sigma_2 > 0.0)
            & (
                phase_II["colorK_36_0"] - sigma_2
                > -2.85714 * (phase_II["color36_45_0"] - sigma_1 - 0.101) + 0.5
            )
        )
    ]

    class_II = safe_setdiff(class_II, phase_II, keys="designation")

    _, ax = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

    ax[0].scatter(
        phase_II["h_m"] - phase_II["k_m"],
        phase_II["j_m"] - phase_II["h_m"],
        s=2,
        color="grey",
    )
    ax[0].scatter(
        class_II["h_m"] - class_II["k_m"],
        class_II["j_m"] - class_II["h_m"],
        s=3,
        color="black",
    )
    ax[0].set_ylabel("J-H")
    ax[0].set_xlabel("H-K")
    ax[0].set_xlim([-0.6, 5])
    ax[0].set_ylim([0, 5])

    ax[1].scatter(phase_II["color36_45_0"], phase_II["colorHK_0"], s=2, color="grey")
    ax[1].scatter(class_II["color36_45_0"], class_II["colorHK_0"], s=3, color="black")
    ax[1].set_ylabel("[H-K]_0")
    ax[1].set_xlabel("[[3.6]-[4.5]]_0")
    ax[1].set_xlim([-1, 1.5])
    ax[1].set_ylim([-0.5, 2])

    plt.savefig("../output/figures/YSO_phase_II_fig17_" + cat + ".svg", format="svg")

    plt.show()

    fig, ax = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)

    ax[0].scatter(phase_II["color36_45_0"], phase_II["colorHK_0"], s=2, color="grey")
    ax[0].scatter(class_II["color36_45_0"], class_II["colorHK_0"], s=5, color="black")

    ax[0].set_ylabel("[H-K]_0")
    ax[0].set_xlabel("[[3.6]-[4.5]]_0")
    # ax[0].set_xlim([-1, 1.5])
    # ax[0].set_ylim([-0.5, 2])

    ax[1].scatter(
        phase_II["color36_45"],
        phase_II["h_m"] - phase_II["k_m"],
        s=2,
        color="grey",
    )
    ax[1].scatter(
        class_II["color36_45"],
        class_II["h_m"] - class_II["k_m"],
        s=5,
        color="black",
    )
    ax[1].set_ylabel("[H-K]")
    ax[1].set_xlabel("[3.6 - 4.5]")
    # ax[1].set_xlim([-0.6, 5])
    # ax[1].set_ylim([0, 5])

    fig.savefig("../output/figures/fig17_comp_" + cat + ".svg", format="svg")

    fig = plt.figure()

    class_I = class_II.copy()

    sigma_1 = class_II["d3_6m"] + class_II["d4_5m"]
    sigma_2 = class_II["k_msigcom"] + class_II["d3_6m"]

    class_II = class_II[
        ~(
            class_II["colorK_36_0"] - sigma_2
            > -2.85714 * (class_II["color36_45_0"] - sigma_1 - 0.401) + 1.7
        )
    ]

    class_I = safe_setdiff(class_I, class_II, keys="designation")

    no_class = safe_setdiff(no_class, class_II, keys="designation")

    no_class = safe_setdiff(no_class, class_I, keys="designation")

    print(
        f"[Phase II] - {len(class_II)} Class II objects, {len(class_I)} Class I objects"
    )

    plt.scatter(class_II["color36_45_0"], class_II["colorK_36_0"], color="blue")

    plt.scatter(class_I["color36_45_0"], class_I["colorK_36_0"], color="red")

    # plt.xlim([-1,1.5])
    # plt.ylim([-0.5,2.5])

    plt.ylabel("[K-[3.6]]_0")
    plt.xlabel("[[3.6]-[4.5]]_0")

    plt.savefig(
        "../output/figures/YSO_phase_II_INCOMPLETE_" + cat + ".svg", format="svg"
    )

    plt.show()
    print("AAAA")

    # Testing

    return class_I, class_II, no_class
