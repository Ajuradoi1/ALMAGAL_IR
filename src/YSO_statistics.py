import pandas as pd
from matplotlib import pyplot as plt
from astropy.table import Table, setdiff,join
from astropy import units as u
import numpy as np
import os


#def YSO_plots(spitzer_sources,alma_regions,alma_sources, alma_proceeding,cat):

def YSO_plots_compare(spitzer_sources,alma_regions,alma_sources,spitzer_sources_alma,alma_regions_alma,alma_sources_alma,alma_proceeding,cat):

    # pah_galaxies|agn_sources|sed_sources|pah_contaminated|p1_class_I|p1_class_II|p1_noclass|
    # p2_class_I|p2_class_II|p2_noclass|

    # MARK: Histogram of sources

    # Spitzer

    cols = ["pah_galaxies",
        "agn_sources",
        "sed_sources",
        "pah_contaminated",
        "p1_class_I",
        "p1_class_II",
        "p1_noclass",
    ]

    fig, axs = plt.subplots(1,2,figsize=(5,4),tight_layout=True,gridspec_kw={'width_ratios': [7, 3]},sharey=True)

    fig.subplots_adjust(bottom=0.25)

    for idx, table in enumerate([spitzer_sources,spitzer_sources_alma]):

        counts = []

        for col in cols:
            
            counts.append(len(table[col][table[col] == 'True']))

        x = np.arange(len(counts))

        if idx == 0:
            bars = axs[0].bar(x, counts,color='lightgrey')
        else:
            bars = axs[0].bar(x, counts,color='magenta')

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
        
        if idx == 0:

            max_phase_1 = np.max(counts)
            axs[0].set_ylim(1,3*max_phase_1)

    axs[0].set_xticks(x)
    axs[0].set_xticklabels(['PAH Gal.','AGN','SED', 'PAH Contam.', 'Class I', 'Class II', 'No Class'], rotation=45, ha="right", fontsize=10)
    axs[0].set_yscale('log')

    

    counts = []

    cols_2 = ["p2_class_I", "p2_class_II", "p2_noclass"]
    for col in cols_2:
        counts.append(len(table[col][table[col] == 'True']))

    x = np.arange(len(counts))

    bars = axs[1].bar(x, counts,color='magenta')

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
    axs[1].set_xticklabels(['Class I', 'Class II', 'No Class'], rotation=45, ha="right", fontsize=10)
    axs[1].tick_params(axis='y', left=False, right=True, labelright=True, labelleft=False)
    axs[1].tick_params(axis='y', which='minor',left=False, right=True, labelright=True, labelleft=False)

    axs[0].set_ylabel('Sources')

    fig.savefig("../output/figures/"+cat+"/YSO_classification_histogram_"+cat+".pdf",format='pdf',dpi=300)


    # MARK: LUM_MASS

    spitzer_sources = spitzer_sources[((spitzer_sources['p1_class_I'] == 'True') | (spitzer_sources['p1_class_II'] == 'True')) ]

    spitzer_sources_alma = spitzer_sources_alma[((spitzer_sources_alma['p1_class_I'] == 'True') | (spitzer_sources_alma['p1_class_II'] == 'True') | (spitzer_sources_alma['p2_class_I'] == 'True') |(spitzer_sources_alma['p2_class_II'] == 'True'))]

    fig, axs = plt.subplots(2,1,figsize=(5,4),tight_layout=True,sharex=True)

    for idx, table in enumerate([spitzer_sources,spitzer_sources_alma]):

        bins_array = np.array(range(-15,33,3)) / 10

        if idx == 0:
            class_I = table['LUM'][(table['p1_class_I'] == 'True')] / table['MASS'][(table['p1_class_I'] == 'True')]
            class_II = table['LUM'][(table['p1_class_II'] == 'True')] / table['MASS'][(table['p1_class_II'] == 'True')]
            #logbins = np.logspace(np.log10(np.min([np.min(class_I),np.min(class_II)])),np.log10((np.max([np.max(class_I),np.max(class_II)]))),bins+1)
            axs[0].hist([np.log10(class_I),np.log10(class_II)],stacked=True,color=['blue','red'],bins=bins_array)
            
            #ax.set_ylim(1,10 ** len(str(np.max(np.histogram(counts,bins = range(0,11,1))[0]))))
        else:
            
            class_I = table['LUM'][(table['p1_class_I'] == 'True') | (table['p2_class_I'] == 'True')] / table['MASS'][(table['p1_class_I'] == 'True') | (table['p2_class_I']== 'True')]
            class_II = table['LUM'][(table['p1_class_II'] == 'True') | (table['p2_class_II']== 'True')] / table['MASS'][(table['p1_class_II'] == 'True') | (table['p2_class_II']== 'True')]
            full_L_M = alma_sources['Lclump'] / alma_sources['Mclump']
            #logbins = np.logspace(np.log10(np.min([np.min(class_I),np.min(class_II)])),np.log10((np.max([np.max(class_I),np.max(class_II)]))),bins+1)
            axs[1].hist(np.log10(full_L_M),stacked=True,color='lightgrey',bins=bins_array)
            axs[1].hist([np.log10(class_I),np.log10(class_II)],stacked=True,color=['blue','red'],bins=bins_array)

        axs[1].set_xlabel('L/M of sources')
        axs[0].set_ylabel('Spitzer sources')
        axs[1].set_ylabel('ALMAGAL matches')

        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        #axs[0].set_xscale('log')
        #axs[1].set_xscale('log')

        fig.savefig('../output/figures/'+cat+'/YSO_'+cat+'_LUM_MASS.pdf',format='pdf',dpi=300)

    
    #MARK: n(H2)

    SURFD_to_nH2 = (
            
            6.022e23  # Avogadro's constant (N_A) [particles/mol]
            / 2.016  # molecular mass of n(H2) [g/mol]
        )

    fig, axs = plt.subplots(2,1,figsize=(4,6),tight_layout=True,sharey=True,sharex=True)

    for idx, table in enumerate([spitzer_sources,spitzer_sources_alma]):

        
        bins_array = np.array(range(2160,2520,20)) / 100
        if idx == 0:
            class_I_DECONV = table['DECONV_SURFD'][(table['p1_class_I'] == 'True')] * SURFD_to_nH2
            class_II_DECONV = table['DECONV_SURFD'][(table['p1_class_II'] == 'True')] * SURFD_to_nH2

            class_I_NODEC = table['NODEC_SURFD'][(table['p1_class_I'] == 'True')] * SURFD_to_nH2
            class_II_NODEC = table['NODEC_SURFD'][(table['p1_class_II'] == 'True')] * SURFD_to_nH2

            #logbins = np.logspace(np.log10(np.min([np.min(class_I),np.min(class_II)])),np.log10((np.max([np.max(class_I),np.max(class_II)]))),bins+1)
            #axs[0].hist([np.log10(class_I_DECONV),np.log10(class_II_DECONV)],stacked=True,color=['blue','red'],bins = bins_array)
            axs[0].hist([np.log10(class_I_NODEC),np.log10(class_II_NODEC)],stacked=True,color=['blue','red'],bins = bins_array)
            
            #ax.set_ylim(1,10 ** len(str(np.max(np.histogram(counts,bins = range(0,11,1))[0]))))
        else:
            class_I_DECONV = table['DECONV_SURFD'][(table['p1_class_I'] == 'True') | (table['p2_class_I'] == 'True')] *SURFD_to_nH2
            class_II_DECONV = table['DECONV_SURFD'][(table['p1_class_II'] == 'True') | (table['p2_class_II']== 'True')] *SURFD_to_nH2

            class_I_NODEC = table['NODEC_SURFD'][(table['p1_class_I'] == 'True') | (table['p2_class_I'] == 'True')] *SURFD_to_nH2
            class_II_NODEC = table['NODEC_SURFD'][(table['p1_class_II'] == 'True') | (table['p2_class_II']== 'True')] *SURFD_to_nH2

            full_NODEC = alma_sources['Surfd']*SURFD_to_nH2
            #logbins = np.logspace(np.log10(np.min([np.min(class_I),np.min(class_II)])),np.log10((np.max([np.max(class_I),np.max(class_II)]))),bins+1)
            #axs[1][0].hist([np.log10(class_I_DECONV),np.log10(class_II_DECONV)],stacked=True,color=['blue','red'],bins = bins_array)
            axs[1].hist([np.log10(full_NODEC)],color='lightgrey',bins=bins_array)
            axs[1].hist([np.log10(class_I_NODEC),np.log10(class_II_NODEC)],stacked=True,color=['blue','red'],bins = bins_array)
            
                
        for i in range(2):
            axs[i].set_yscale('log')
            #axs[i].set_ylim(0,1e4)
        
        axs[0].set_ylabel('Spitzer sources')
        axs[1].set_ylabel('ALMAGAL matches')

        fig.savefig('../output/figures/'+cat+'/YSO_'+cat+'_NH2.pdf',format='pdf',dpi=300)

def region_plot(spitzer_sources,spitzer_sources_alma,cat,region):

    fig, ax = plt.subplots(figsize=(6,6),tight_layout=True)

    ax.set_xlim([-0.03,0.03])
    ax.set_ylim([-0.03,0.03])

    center = [float(spitzer_sources['l'][spitzer_sources['internal_id_01'] == region][0])*u.deg,float(spitzer_sources['b'][spitzer_sources['internal_id_01'] == region][0])*u.deg]

    for source in spitzer_sources[spitzer_sources['internal_id_01'] == region]:
        ax.scatter(source['l']*u.deg-center[0],source['b']*u.deg-center[1],color='grey')



    os.makedirs('../output/figures/'+cat+'/regions',exist_ok=True)
    fig.savefig('../output/figures/'+cat+'/regions/'+region+'.pdf',format='pdf',dpi=300)

    return 0








    



