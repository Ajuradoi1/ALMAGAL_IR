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

from spitzer_sources_script import spitzer_sources

from regions_ds9_script import regions_ds9
from clumps_analysis_script import (
    match_table_fun,
    plot_separation_histogram,
    plot_source_matches_all,
    plot_match_heatmap,
)

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


# Loading scripts


if not (
    os.path.isfile("../output/tables/alma_regions_coordinates.tbl")
    and os.path.isfile("../output/tables/spitzer_sources_" + sys.argv[1] + ".tbl")
):
    print("\ndownloading spitzer catalog...\n")
    spitzer_sources(sys.argv[2], sys.argv[1])
else:
    print("\nspitzer catalog found, skipping download...")

print("\ngenerating regions files...\n")

regions_ds9(sys.argv[2], sys.argv[3], sys.argv[1])

print("\nfinding matches...\n")

if os.path.isfile("../output/tables/match_table_" + sys.argv[1] + ".tbl"):
    match_table = Table.read(
        "../output/tables/match_table_" + sys.argv[1] + ".tbl", format="ipac"
    )
else:
    match_table_fun(sys.argv[2], sys.argv[3], sys.argv[1])
    match_table = Table.read(
        "../output/tables/match_table_" + sys.argv[1] + ".tbl", format="ipac"
    )

print("\nplotting...\n")

plot_separation_histogram(match_table, 20, sys.argv[1])

plot_source_matches_all(match_table, sys.argv[1])

plot_match_heatmap(match_table, sys.argv[1])


print("\ndone!")
