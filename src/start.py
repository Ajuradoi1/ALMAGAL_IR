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
from astropy.table import Table

from spitzer_data import spitzer_data

from regions_ds9_script import regions_ds9
from clumps_analysis import (
    match_table_fun,
    plot_separation_histogram,
    plot_source_matches_all,
    plot_match_heatmap,
)

from YSO_Gutermuth_09 import YSO_Phase_I

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
alma_regions = Table.read("../alma_input/" + alma_file, format="ascii.ipac")[:2]

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
    match_table_fun(
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
    print(
        f"\nmatch table found, skipping matching... {len(match_table)} matches found\n"
    )


if (
    input("Do you want to perform the plots from the match table? (y/n) ").lower()
    == "y"
):
    print("\nplotting...\n")

    plot_separation_histogram(match_table=match_table, bins=20, catalog=cat)

    plot_source_matches_all(match_table=match_table, catalog=cat)

    plot_match_heatmap(match_table=match_table, catalog=cat)

    print("\ndone!")

if input("Do you want to perform a YSO classification? (y/n) ").lower() == "y":
    class_I, class_II, no_class = YSO_Phase_I(cat)

    if not os.path.isfile("../output/tables/spitzer_2MASS_fp_psc.tbl"):
        print("\n2MASS data missing, downloading...")
        spitzer_data(alma_regions_file=alma_file, catalog="fp_psc", type="2MASS")

    # spitzer_sources = spitzer_sources[spitzer_sources["designation"].isin(match_table["spitzer_source"])]
