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
from astropy.table import join
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

from YSO_statistics import YSO_plots_compare

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


# Loading scripts


if not os.path.isfile("../output/tables/"+cat+"/spitzer_sources_" + cat + "_download.tbl"):
    print("\ndownloading spitzer catalog...\n")
    spitzer_data(alma_regions_file=alma_file, catalog=cat, type="sources")
    spitzer_bool = True
else:
    print("\nthis catalog already exists, skipping download...")
    spitzer_bool = False


spitzer_sources = Table.read(
    "../output/tables/"+cat+"/spitzer_sources_" + cat + "_download.tbl",
    format="ipac",
)
spitzer_sources = unique(spitzer_sources, keys="designation")

obs_sample_map = dict(zip(alma_regions["INTERNAL_ID"], alma_regions["OBS_SAMPLE"]))

obs_array = np.array(
    [obs_sample_map.get(region, np.nan) for region in spitzer_sources["internal_id_01"]]
)

spitzer_sources["obs_sample"] = obs_array

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
        spitzer_sources=spitzer_sources,alma_regions=alma_regions, alma_sources=alma_sources, catalog=cat
    )
    match_table = Table.read(
        "../output/tables/match_table_" + cat + ".tbl", format="ipac"
    )
    print(f"\n match table generated, {len(match_table)} matches found\n")
    match_bool = True
else:
    match_table = Table.read(
        "../output/tables/match_table_" + cat + ".tbl", format="ipac"
    )
    print(f"\nmatch table found, skipping matching... {len(match_table)} matches found")
    print(f"\n{len(set(match_table['region_id']))} ALMA regions have matches\n")
    match_bool = False

spitzer_sources["has_alma"] = [
    i in (match_table["spitzer_source"]) for i in spitzer_sources["designation"]
]

#if spitzer_bool:
if True:
    spitzer_sources = join(
                spitzer_sources,
                match_table["region_id", "alma_source", "spitzer_source"],
                keys_left = ['designation'],
                keys_right = ["spitzer_source"],
                join_type = "left",
            )

    spitzer_sources.remove_columns(["region_id","spitzer_source"])

alma_regions["HAS_SPITZER"] = [
    r in match_table['region_id'] for r in alma_regions['INTERNAL_ID']
]

'''alma_regions = join(
            alma_regions,
            match_table["region_id", "alma_source", "spitzer_source"],
            keys_left = ["INTERNAL_ID"],
            keys_right = ["region_id"],
            join_type = "left",
        )

alma_regions.remove_columns(["region_id"])'''

match_keys_alma = set(zip(match_table["region_id"], match_table["alma_source"]))

alma_sources["HAS_SPITZER"] = [
    (r, s) in match_keys_alma for r, s in zip(alma_sources["CLUMP"], alma_sources["ID"])
]

alma_sources = join(
            alma_sources,
            match_table["region_id", "alma_source", "spitzer_source"],
            keys_left = ["CLUMP", "ID"],
            keys_right = ["region_id", "alma_source"],
            join_type = "left",
        )

alma_sources.remove_columns(["region_id", "alma_source"])


if match_bool:

    print('\ngenerating '+cat+' with ALMA matches table...\n')

    spitzer_sources.write(
        "../output/tables/"+cat+"/spitzer_sources_" + cat + ".tbl", format="ascii.ipac", overwrite=True
    )

    print('\ngenerating alma regions with '+cat+' matches table...\n')

    alma_regions.write(
        "../output/tables/"+cat+"/ALMA_REGIONS_" + cat + ".tbl", format="ascii.ipac", overwrite=True
    )

    print('\ngenerating alma sources with '+cat+' matches table...\n')

    alma_sources.write(
        "../output/tables/"+cat+"/ALMA_SOURCES_" + cat + ".tbl", format="ascii.ipac", overwrite=True
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
    cat_alma = cat + "_ALMA"
    alma_proceeding = True
else:
    cat_alma = cat
    alma_proceeding = False

if input("Do you want to perform a YSO classification? (y/n) ").lower() == "y":
    # We need to add 2MASS data to our spitzer catalog

    if not os.path.isfile("../output/tables/spitzer_2MASS_fp_psc_download.tbl"):
        print("\n2MASS data missing, downloading...")
        spitzer_data(alma_regions_file=alma_file, catalog="fp_psc", type="2MASS")

    print("\nadding 2mass data to spitzer...")

    spitzer_sources = twomass_to_spitzer(spitzer_sources, cat)
    

    cols_spitzer = ['Ak_DECONV', 'Ak_NODEC']

    for col in cols_spitzer:
        spitzer_sources[col] = 0

    spitzer_sources = join(
        spitzer_sources,
        alma_regions['INTERNAL_ID','MASS','LUM','DECONV_SURFD','NODEC_SURFD'],
        keys_left = ['internal_id_01'],
        keys_right = ['INTERNAL_ID'],
        join_type = "left",
    )
      
    spitzer_sources.remove_columns(['INTERNAL_ID'])

    SURFD_to_Ak = (
            (
                0.112 / 2.21e21
            )  # Ak = 0.112 Av, Av = (2.21e21)^-1 n(H2). C_k = 5.07e-23 [mag cm^2]
            * 6.022e23  # Avogadro's constant (N_A) [particles/mol]
            / 2.016  # molecular mass of n(H2) [g/mol]
        )


    if alma_proceeding:
        

        spitzer_sources = join(
            spitzer_sources,
            alma_sources['spitzer_source','Mclump','Lclump','Lclump_Mclump'][alma_sources['HAS_SPITZER']],
            keys_left = ['designation'],
            keys_right = ['spitzer_source'],
            join_type = "left",
        )
        
        spitzer_sources.remove_columns(["spitzer_source"])


    # Add columns
    spitzer_sources["Ak_DECONV"] = spitzer_sources['DECONV_SURFD'] * SURFD_to_Ak
    spitzer_sources["Ak_NODEC"] = spitzer_sources['NODEC_SURFD'] * SURFD_to_Ak
    mean_ak_nodec = np.mean(
        spitzer_sources["Ak_NODEC"][
            spitzer_sources["Ak_NODEC"] != 0
        ]
    )
    mean_ak_deconv = np.mean(
        spitzer_sources["Ak_DECONV"][
            spitzer_sources["Ak_DECONV"] != 0
        ]
    )

    for source in spitzer_sources:
        if source["Ak_NODEC"] == 0:
            source["Ak_NODEC"] = mean_ak_nodec
        if source["Ak_DECONV"] == 0:
            source["Ak_DECONV"] = mean_ak_deconv

    # Plot Ak hist

    fig = plt.figure(figsize=(4, 5), tight_layout=True)
    fig.subplots_adjust(bottom=0.25)
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
    plt.savefig("../output/figures/"+cat+"/Ak_histogram_" + cat_alma + ".pdf", format="pdf")

    spitzer_sources.write(
    "../output/tables/"+cat+"/spitzer_sources_" + cat_alma + "_2mass.tbl",
    format="ascii.ipac",
    overwrite=True,
    )  

    spitzer_sources = YSO_Phase_I(spitzer_sources, cat)

    

    if alma_proceeding:

        spitzer_sources.write(
        "../output/tables/"+cat+"/YSO_phase_I_" + cat_alma + ".tbl",
        format="ascii.ipac",
        overwrite=True,
        )


        print("\nproceeding with phase II...")
        spitzer_sources = YSO_Phase_II(spitzer_sources, cat)

        spitzer_sources.write(
            "../output/tables/"+cat+"/YSO_phase_II_" + cat_alma + ".tbl",
            format="ascii.ipac",
            overwrite=True,
        )
        print("\nphase II completed")


        alma_sources = join(alma_sources,
                            spitzer_sources["internal_id_01",
                                            "alma_source",
                                            "designation",
                                            "pah_galaxies",
                                            "agn_sources",
                                            "sed_sources",
                                            "pah_contaminated",
                                            "p1_class_I",
                                            "p1_class_II",
                                            "p1_noclass",
                                            "p2_class_I",
                                            "p2_class_II",
                                            "p2_noclass"],
                            keys_left = ['CLUMP','ID'],
                            keys_right= ['internal_id_01','alma_source'],
                            join_type='left'
                            )


        

        alma_sources.remove_columns(['internal_id_01','alma_source'])


        spitzer_sources.write('../output/tables/'+cat+'/spitzer_sources_YSO_PLOTS_ALMA.tbl', format='ascii.ipac', overwrite=True)

        alma_sources.write('../output/tables/'+cat+'/alma_sources_YSO_PLOTS_ALMA.tbl', format='ascii.ipac', overwrite=True)

        alma_regions.write('../output/tables/'+cat+'/alma_regions_YSO_PLOTS_ALMA.tbl', format='ascii.ipac', overwrite=True)

    else:

        spitzer_sources.write(
        "../output/tables/"+cat+"/YSO_phase_I_" + cat + ".tbl",
        format="ascii.ipac",
        overwrite=True,
        )

        print("\nphase II is not implemented for sources without ALMA matches")

        alma_sources_join = join(alma_sources['CLUMP','ID'][alma_sources['HAS_SPITZER']],
                            spitzer_sources["internal_id_01",
                                            "alma_source",
                                            "designation",
                                            "pah_galaxies",
                                            "agn_sources",
                                            "sed_sources",
                                            "pah_contaminated",
                                            "p1_class_I",
                                            "p1_class_II",
                                            "p1_noclass"][spitzer_sources['has_alma']],
                            keys_left = ['CLUMP','ID'],
                            keys_right= ['internal_id_01','alma_source'],
                            join_type='left'
                            )
        alma_sources_join.rename_column('CLUMP','CLUMP_join')
        alma_sources_join.rename_column('ID','ID_join')

        alma_sources = join(
            alma_sources,
            alma_sources_join,
            keys_left = ['CLUMP','ID'],
            keys_right = ['CLUMP_join','ID_join'],
            join_type='left'
        )

        alma_sources.remove_columns(['CLUMP_join','ID_join'])

        alma_sources.remove_columns(['internal_id_01','alma_source'])

        spitzer_sources.write('../output/tables/'+cat+'/spitzer_sources_YSO_PLOTS.tbl', format='ascii.ipac', overwrite=True)

        alma_sources.write('../output/tables/'+cat+'/alma_sources_YSO_PLOTS.tbl', format='ascii.ipac', overwrite=True)

        alma_regions.write('../output/tables/'+cat+'/alma_regions_YSO_PLOTS.tbl', format='ascii.ipac', overwrite=True)
else:
    print("\nanalysis finished.")


# YSO_plots(spitzer_sources, alma_proceeding, cat)
