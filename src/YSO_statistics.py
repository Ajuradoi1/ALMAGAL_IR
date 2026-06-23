import pandas as pd
from matplotlib import pyplot as plt
from astropy.table import Table, setdiff,join
import numpy as np
import os


def YSO_plots(spitzer_sources,alma_regions,alma_sources, alma_proceeding,cat):

    # pah_galaxies|agn_sources|sed_sources|pah_contaminated|p1_class_I|p1_class_II|p1_noclass|
    # p2_class_I|p2_class_II|p2_noclass|

    # MARK: COUNT

    cols = ["designation",
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
        alma_sources = alma_sources[alma_sources['HAS_SPITZER'] == 'True'] #change in real code

        alma_sources = join(
            alma_sources,
            spitzer_sources[cols],
            keys_left = ['spitzer_source'],
            keys_right = ['designation'],
            join_type = "left",
        )

        cols.pop(0)

        alma_sources.remove_columns(['designation'])

        counts = []

        print(cols)

        for col in cols:
            counts.append(len(alma_sources[col][alma_sources[col] == 'True']))
            print(f"\nNumber of {col}: {counts[-1]}\n")

        '''idx = np.argsort(counts)[::-1]

        cols = np.array(cols)[idx]
        counts = np.array(counts)[idx]'''


        fig, axs = plt.subplots(1,2,figsize=(5,4),tight_layout=True,gridspec_kw={'width_ratios': [7, 3]},sharey=True)
        fig.subplots_adjust(bottom=0.25)
        print(len(counts))
        x = np.arange(len(counts[:7]))
        bars = axs[0].bar(x, counts[:7])

        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            axs[0].text(
                bar.get_x() + bar.get_width() / 2,
                height,
                f"{height}",
                ha="center",
                va="bottom",
                fontsize=9
            )

        axs[0].set_xticks(x)
        axs[0].set_xticklabels(['PAH Gal.','AGN','SED', 'PAH Contam.', 'Class I', 'Class II', 'Unclassified'], rotation=45, ha="right", fontsize=10)

        max_phase_1 = sum(counts[:7])
        max_phase_2 = sum(counts[7:])

        x = np.arange(len(counts[7:]))
        bars = axs[1].bar(x, counts[7:])

        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            axs[1].text(
                bar.get_x() + bar.get_width() / 2,
                height,
                f"{height}",
                ha="center",
                va="bottom",
                fontsize=9
            )

        axs[1].set_xticks(x)
        axs[1].set_xticklabels(['Class I', 'Class II', 'Unclassified'], rotation=45, ha="right", fontsize=10)
        axs[1].tick_params(axis='y', left=False, right=True, labelright=True, labelleft=False)

        

        #plt.legend()
        
        axs[0].set_ylim(1,2*max_phase_1)
        #axs[1].set_ylim(1,max_phase_2)

        axs[0].set_yscale('log')
        axs[1].set_yscale('log')

        fig.savefig("../output/figures/YSO_classification_histogram_ALMA.pdf",format='pdf',dpi=300)

    else:

        cols.pop(0)

        fig, ax = plt.subplots(figsize=(5,4),tight_layout=True)
        fig.subplots_adjust(bottom=0.25)
        print(cols)

        counts = []

        for col in cols:
            counts.append(len(spitzer_sources[col][spitzer_sources[col] == 'True']))
            print(f"\nNumber of {col}: {counts[-1]}\n")


        print(len(counts))
        

        x = np.arange(len(counts))
        bars = ax.bar(x, counts)

        # Add value labels on top of bars
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height,
                f"{height}",
                ha="center",
                va="bottom",
                fontsize=9
            )

        ax.set_xticks(x)
        ax.set_xticklabels(['PAH Gal.','AGN','SED', 'PAH Contam.', 'Class I', 'Class II', 'Unclassified'], rotation=45, ha="right", fontsize=10)
        ax.set_yscale('log')
        max_phase_1 = sum(counts)
        ax.set_ylim(0,2*max_phase_1)

        fig.savefig("../output/figures/YSO_classification_histogram_glimpse_s07.pdf",format='pdf',dpi=300)

