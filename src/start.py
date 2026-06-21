# ALMAGAL Source analyzer -- Spitzer analyzer
# Author: Álex Jurado i López
# github repo: https://github.com/Ajuradoi1/TFM-ALMA

# Execute the script giving <spitzer_catalog> <alma_regions_file>(opt) <alma_source_file>(opt)

# For testing purposes:
# <alma_regions_file> almaprops_2024_all_20250219.tab
# <alma_source_file> 7MTM2TM1_Core_catalogue_out_official_v3_run6Apr23+16Jan24_sn-5_may24_clump_cat.txt
# <spitzer_catalog> glimpse_s07 // glimpse2_v2cat

import sys
import os
from astropy.table import Table, unique
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd
from spitzer_data import spitzer_data

from regions_ds9_script import regions_ds9
from clumps_analysis import (
    alma_match_table_fun,
    plot_separation_histogram,
    plot_source_matches_all,
    plot_match_heatmap,
)

from twomass_to_spitzer import twomass_to_spitzer

from YSO_Gutermuth_09 import YSO_Phase_I, YSO_Phase_II

from YSO_statistics import YSO_plots

if len(sys.argv) == 1:
    sys.argv.append("glimpse_s07")
    sys.argv.append("almaprops_2024_all_20250219.tab")
    sys.argv.append(
        "7MTM2TM1_Core_catalogue_out_official_v3_run6Apr23+16Jan24_sn-5_may24_clump_cat.txt"
    )
elif len(sys.argv) == 2:
    sys.argv.append("almaprops_2024_all_20250219.tab")
    sys.argv.append(
        "7MTM2TM1_Core_catalogue_out_official_v3_run6Apr23+16Jan24_sn-5_may24_clump_cat.txt"
    )

cat = sys.argv[1]
alma_file = sys.argv[2]
alma_sources_file = sys.argv[3]
alma_regions = Table.read("../alma_input/" + alma_file, format="ascii.ipac")


# alma_sources = Table.read("../alma_input/" + alma_sources_file, format="ascii.ipac")
alma_sources = Table.from_pandas(
    pd.read_csv(
        "../alma_input/" + alma_sources_file,
        comment="\\",
        skiprows=[2],
        delimiter=r"\s+",
    )
)

alma_sources.rename_column("S/N", "S_N")
alma_sources.rename_column("FWHM_max/FWHM_min", "FWHM_max__FWHM_min")
alma_sources.rename_column("Bmaj/Bmin", "Bmaj_Bmin")

alma_sources.rename_column("CLUMP.1", "CLUMP_1")

alma_sources.rename_column("Dcore/BEAM_size", "Dcore__Beam_size")
alma_sources.rename_column("Lclump/Mclump", "Lclump_Mclump")
alma_sources.rename_column("n(H2)", "n_H2")


# alma_regions = alma_regions[0]

# Loading scripts


if not os.path.isfile("../output/tables/spitzer_sources_" + cat + ".tbl"):
    print("\ndownloading spitzer catalog...\n")
    spitzer_data(alma_regions_file=alma_file, catalog=cat, type="sources")
else:
    print("\nthis catalog already exists, skipping download...")


if (
    input("\ndoes this catalog has archival data you want to add? (y/n) ").lower()
    == "y"
):
    arch = input("\nplease enter the archival data catalog:")
    cols = (
        input('\nwhat columns do you want to add? (enter them separated by ",")')
        .lower()
        .split(",")
    )

    spitzer_table = Table.read(
        "../output/tables/spitzer_sources_" + cat + ".tbl",
        format="ipac",
    )
    arch_table = Table.read("../output/tables/spitzer_sources_" + arch + ".tbl")

    for col in cols:
        spitzer_table = spitzer_table.add_column(
            arch_table[spitzer_table["designation"] == arch_table["designation"]][col],
            name=col,
            index=-1,
        )
    spitzer_table = spitzer_table.drop_duplicates(subset="designation")
    spitzer_table.write(
        "../output/tables/spitzer_sources_" + cat + ".tbl",
        format="ascii.ipac",
        overwrite=True,
    )

    print(f"\n your table is ready! {len(cols)} columns have been added")
else:
    print("\nskipping...")

if not os.path.isfile(
    "../output/alma_sources/"
    + str(alma_regions["INTERNAL_ID"][0])
    + "/"
    + str(alma_regions["INTERNAL_ID"][0])
    + "_"
    + cat
    + "_spitzer_sources.reg"
):
    print("\ngenerating regions files...\n")
    regions_ds9(alma_file, alma_sources_file, cat)
else:
    print("\nregions file found...")

print("\nfinding matches...\n")

if not os.path.isfile("../output/tables/match_table_" + cat + ".tbl"):
    alma_match_table_fun(
        alma_regions_file=alma_file, alma_sources_file=alma_sources_file, catalog=cat
    )
    match_table = Table.read(
        "../output/tables/match_table_" + cat + ".tbl", format="ipac"
    )
    print(f"\n match table generated, {len(match_table)} matches found\n")
else:
    match_table = Table.read(
        "../output/tables/match_table_" + cat + ".tbl", format="ipac"
    )
    print(f"\nmatch table found, skipping matching... {len(match_table)} matches found")
    print(f"\n{len(set(match_table['region_id']))} ALMA regions have matches\n")

spitzer_sources = Table.read(
    "../output/tables/spitzer_sources_" + cat + ".tbl",
    format="ipac",
)
spitzer_sources = unique(spitzer_sources, keys="designation")

spitzer_sources["has_alma"] = [
    i in (match_table["spitzer_source"]) for i in spitzer_sources["designation"]
]

match_keys_alma = set(zip(match_table["region_id"], match_table["alma_source"]))

alma_sources["HAS_SPITZER"] = [
    (r, s) in match_keys_alma for r, s in zip(alma_sources["CLUMP"], alma_sources["ID"])
]


alma_sources.write(
    "../output/tables/ALMA_SOURCES_" + cat + ".tbl", format="ascii.ipac", overwrite=True
)

if (
    input("Do you want to perform the plots from the match table? (y/n) ").lower()
    == "y"
):
    print("\nplotting...\n")

    plot_separation_histogram(match_table=match_table, bin_width=5, catalog=cat)

    plot_source_matches_all(match_table=match_table, catalog=cat)

    plot_match_heatmap(match_table=match_table, catalog=cat)

    print("\ndone!")

if (
    input("Do you want to work onwards with ALMA matched sources? (y/n) ").lower()
    == "y"
):
    spitzer_sources = spitzer_sources[(spitzer_sources["has_alma"])]
    cat_name = cat
    cat = cat + "_ALMA"
    alma_proceeding = True
else:
    cat_name = cat
    alma_proceeding = False

if input("Do you want to perform a YSO classification? (y/n) ").lower() == "y":
    # We need to add 2MASS data to our spitzer catalog

    if not os.path.isfile("../output/tables/spitzer_sources_" + cat + "_2mass.tbl"):
        if not os.path.isfile("../output/tables/spitzer_2MASS_fp_psc.tbl"):
            print("\n2MASS data missing, downloading...")
            spitzer_data(alma_regions_file=alma_file, catalog="fp_psc", type="2MASS")

        print("\nadding 2mass data to spitzer...")

        spitzer_sources = twomass_to_spitzer(spitzer_sources, cat)
    else:
        spitzer_sources = Table.read(
            "../output/tables/spitzer_sources_" + cat + "_2mass.tbl",
            format="ascii.ipac",
        )

    if alma_proceeding:
        # SURFD_to_Ak = (1 / (2.8)) * (6.022 * 1e23) * (1 / (8.4 * 1e21))
        # 25
        # SURFD_to_Ak = 1 / (2.8 * 1.673 * 1e-24 * 8.4 * 1e21)
        SURFD_to_Ak = (
            (
                0.112 / 2.21e21
            )  # Ak = 0.112 Av, Av = (2.21e21)^-1 n(H2). C_k = 5.07e-23 [mag cm^2]
            * 6.022e23  # Avogadro's constant (N_A) [particles/mol]
            / 2.016  # molecular mass of n(H2) [g/mol]
        )
        # Build fast lookups
        deconv_lookup = dict(
            zip(alma_regions["INTERNAL_ID"], alma_regions["DECONV_SURFD"])
        )

        nodec_lookup = dict(
            zip(alma_regions["INTERNAL_ID"], alma_regions["NODEC_SURFD"])
        )

        m_lookup = dict(
            zip(
                zip(alma_sources["CLUMP"], alma_sources["ID"]),
                alma_sources["Mclump"],
            )
        )
        l_lookup = dict(
            zip(
                zip(alma_sources["CLUMP"], alma_sources["ID"]),
                alma_sources["Lclump"],
            )
        )

        mask_m = np.array(
            [
                m_lookup.get((region, source_id), np.nan)
                for region, source_id in zip(
                    match_table["region_id"], match_table["alma_source"]
                )
            ]
        )

        mask_l = np.array(
            [
                l_lookup.get((region, source_id), np.nan)
                for region, source_id in zip(
                    match_table["region_id"], match_table["alma_source"]
                )
            ]
        )

        spitzer_sources["Mclump"] = mask_m
        spitzer_sources["Lclump"] = mask_l

        # Fill arrays safely
        mask_DECONV = np.array(
            [
                deconv_lookup.get(region, np.nan)
                for region in spitzer_sources["internal_id_01"]
            ]
        )

        mask_NODEC = np.array(
            [
                nodec_lookup.get(region, np.nan)
                for region in spitzer_sources["internal_id_01"]
            ]
        )

        # Add columns
        spitzer_sources["Ak_DECONV"] = mask_DECONV * SURFD_to_Ak
        spitzer_sources["Ak_NODEC"] = mask_NODEC * SURFD_to_Ak
        mean_ak_nodec = np.mean(
            spitzer_sources["Ak_NODEC"][
                spitzer_sources["Ak_NODEC"] != -999 * SURFD_to_Ak
            ]
        )
        mean_ak_deconv = np.mean(
            spitzer_sources["Ak_DECONV"][
                spitzer_sources["Ak_DECONV"] != -999 * SURFD_to_Ak
            ]
        )

        for source in spitzer_sources:
            if source["Ak_NODEC"] == -999 * SURFD_to_Ak:
                source["Ak_NODEC"] = mean_ak_nodec
            if source["Ak_DECONV"] == -999 * SURFD_to_Ak:
                source["Ak_DECONV"] = mean_ak_deconv

        fig = plt.figure(figsize=(4, 6), tight_layout=True)

        plt.hist(
            spitzer_sources["Ak_NODEC"],
            bins=range(0, 125, 5),
            alpha=0.5,
            label="NODEC",
            color="blue",
        )
        plt.hist(
            spitzer_sources["Ak_DECONV"],
            bins=range(0, 125, 5),
            alpha=0.5,
            label="DECONV",
            color="red",
        )
        plt.ylabel("Regions")
        plt.xlabel("Ak")
        plt.xlim([0, 125])
        plt.legend()
        plt.savefig("../output/figures/Ak_histogram_" + cat_name + ".pdf", format="pdf")

    spitzer_sources.write(
        "../output/tables/spitzer_sources_" + cat + "_2mass.tbl",
        format="ascii.ipac",
        overwrite=True,
    )

    spitzer_sources = YSO_Phase_I(spitzer_sources, cat)

    spitzer_sources.write(
        "../output/tables/YSO_phase_I_" + cat + ".tbl",
        format="ascii.ipac",
        overwrite=True,
    )

    if alma_proceeding:
        print("\nproceeding with phase II...")
        spitzer_sources = YSO_Phase_II(spitzer_sources, cat)

        spitzer_sources.write(
            "../output/tables/YSO_phase_II_" + cat + ".tbl",
            format="ascii.ipac",
            overwrite=True,
        )
        print("\nphase II completed")

    else:
        print("\nphase II is not implemented for sources without ALMA matches")

else:
    print("\nanalysis finished.")


# YSO_plots(spitzer_sources, alma_proceeding, cat)
