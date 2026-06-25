import pandas as pd
from matplotlib import pyplot as plt
from astropy.table import Table, setdiff,join
from astropy import units as u
import numpy as np
import os
import shapely
from shapely import MultiPoint, LineString, Point
from astropy.coordinates import SkyCoord


#def YSO_plots(spitzer_sources,alma_regions,alma_sources, alma_proceeding,cat):

def YSO_plots_compare(spitzer_sources,alma_sources,spitzer_sources_alma,cat):

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
        axs[1].set_xlabel(r'log $n(H_2)$')

        fig.savefig('../output/figures/'+cat+'/YSO_'+cat+'_NH2.pdf',format='pdf',dpi=300)

#MARK: Region analyzer

def region_analyze(spitzer_sources,alma_regions,cat,alma_proceeding):
    
    #pah_galaxies, agn_sources, sed_sources, pah_contam, p1_class_I, p1_class_II, p1_noclass, p2_class_I, p2_class_II, p2_noclass 
    if alma_proceeding:
        cols = ['pah_galaxies', 'agn_sources', 'sed_sources', 'pah_contaminated','p1_class_I', 'p1_class_II', 'p1_noclass', 'p2_class_I', 'p2_class_II', 'p2_noclass'] 
    else:
        cols = ['pah_galaxies', 'agn_sources', 'sed_sources', 'pah_contaminated','p1_class_I', 'p1_class_II', 'p1_noclass']

    for col in cols:
        alma_regions[col] = 0
    
    for region in alma_regions:
        region_id = region['INTERNAL_ID']

        for col in cols:
            alma_regions[col][alma_regions['INTERNAL_ID'] == region_id] = int(len(spitzer_sources[col][(spitzer_sources[col] == 'True') & (spitzer_sources['internal_id_01'] == region_id)]))
    alma_regions.write('../output/tables/'+cat+'/alma_regions_YSO.tbl',format='ascii.ipac',overwrite=True)

    # region max sources

    max_count_total = 0
    max_count_YSO = 0

    count_dict = {}



    for idx,region in enumerate(alma_regions):
        
        count_total = 0
        count_YSO = 0
        for col in cols:
            count_total = count_total + region[col]
            if alma_proceeding:
                if col in ['p1_class_I','p1_class_II','p2_class_I','p2_class_II']:
                    count_YSO = count_YSO + region[col]
            else:
                 if col in ['p1_class_I','p1_class_II']:
                    count_YSO = count_YSO + region[col]
        
        count_dict[idx] = [region['INTERNAL_ID'],count_total,count_YSO]
    

    total_tol = sorted(count_dict.items(),key=lambda item: item[1][1],reverse=True)
    total_YSO = sorted(count_dict.items(),key=lambda item: item[1][2],reverse=True)

    total_id = []
    yso_id = []

    for element in total_tol:
        total_id.append(element[1][0])
    for element in total_YSO:
        yso_id.append(element[1][0])


    print(f'Region with most sources: {total_id[0]}\n')
    print(f'Region with most YSO matched: {yso_id[0]}\n')

    print(f'\n\n{total_id}\n\n {yso_id}\n\n')
    
    return alma_regions, total_id, yso_id


#MARK: region_plot

def region_plot(spitzer_sources,alma_sources,alma_sources_alma,cat,region,region_info,alma_proceeding,idx):

    fig, axs = plt.subplots(1,2,figsize=(6,6),tight_layout=True)
    axs[0].set_box_aspect(1)
    axs[1].set_box_aspect(1)

    radius = 1 * u.arcmin

    axs[0].set_xlim([radius.to(u.deg).value, -radius.to(u.deg).value])

    axs[0].set_ylim([-radius.to(u.deg).value, radius.to(u.deg).value])

    axs[0].set_xticks(np.linspace(radius.to(u.deg).value,-radius.to(u.deg).value,5))
    axs[0].set_xticklabels(np.linspace(60,-60,5,dtype=int))

    axs[0].set_yticks(np.linspace(-radius.to(u.deg).value,radius.to(u.deg).value,5))
    axs[0].set_yticklabels(np.linspace(-60,60,5,dtype=int))

    axs[0].set_xlabel(r'$\Delta l$ [arcsec]')
    axs[0].set_ylabel(r'$\Delta b$ [arcsec]')

    radius = 0.38/2 * u.arcmin

    axs[1].set_xlim([radius.to(u.deg).value, -radius.to(u.deg).value])

    axs[1].set_ylim([-radius.to(u.deg).value, radius.to(u.deg).value])

    axs[1].set_xticks(np.linspace(radius.to(u.deg).value,-radius.to(u.deg).value,5))
    axs[1].set_xticklabels(np.linspace(19,-19,5,dtype=int))

    axs[1].set_yticks(np.linspace(-radius.to(u.deg).value,radius.to(u.deg).value,5))
    axs[1].set_yticklabels(np.linspace(-19,19,5,dtype=int))

    axs[1].set_xlabel(r'$\Delta l$ [arcsec]')
    axs[1].set_ylabel(r'$\Delta b$ [arcsec]')

    center = [float(region_info['ra'])*u.deg,float(region_info['dec'])*u.deg]

    center_aux = SkyCoord(
        ra =center[0],
        dec = center[1],
        frame='icrs'
    )

    center = [center_aux.galactic.l, center_aux.galactic.b]
    
    # MARK: ALMA_IF

    if alma_proceeding: 

        counter = {}
        
        for i in range(7):
            counter[i] = 0

        class_I = []
        class_II = []

        for i in range(2):
            
            if i == 1:
                class_I = []
                class_II = []

            for idx2, table in enumerate([spitzer_sources,alma_sources,alma_sources_alma]):   
                if idx2 == 0:
                    for source in table:
                        if source['p1_class_I'] == 'True':
                            axs[i].scatter(source['l']*u.deg-center[0],source['b']*u.deg-center[1],color='blue',marker='x',s=50,zorder=4)
                            counter[0] = counter[0] +1
                            class_I.append(((source['l']*u.deg-center[0]).value,(source['b']*u.deg-center[1]).value))
                        elif source['p1_class_II'] == 'True':
                            axs[i].scatter(source['l']*u.deg-center[0],source['b']*u.deg-center[1],color='red',marker='x',s=50,zorder=3)
                            counter[1] = counter[1] +1
                            class_II.append(((source['l']*u.deg-center[0]).value,(source['b']*u.deg-center[1]).value))
                        else:
                            axs[i].scatter(source['l']*u.deg-center[0],source['b']*u.deg-center[1],color='lightgrey',marker='x',s=50,zorder=2)
                            counter[2] = counter[2] +1 
                elif idx2 == 1:
                    for source in table:
                        axs[i].scatter(table['GLON']*u.deg-center[0],table['GLAT']*u.deg-center[1],color='lightgrey',marker='o',s=20,zorder=0)
                        counter[6] = counter[6] +1
                else:
                    for source in table:
                        if ((source['p1_class_I'] == 'True') or (source['p2_class_I'] == 'True')):
                            axs[i].scatter(source['GLON']*u.deg-center[0],source['GLAT']*u.deg-center[1],color='blue',marker='+',s=50,zorder=6)
                            counter[3] = counter[3] +1
                        elif ((source['p1_class_II'] == 'True') or (source['p2_class_II'] == 'True')):
                            axs[i].scatter(source['GLON']*u.deg-center[0],source['GLAT']*u.deg-center[1],color='red',marker='+',s=50,zorder=5)
                            counter[4] = counter[4] +1
                        else:
                            axs[i].scatter(source['GLON']*u.deg-center[0],source['GLAT']*u.deg-center[1],color='lightgrey',marker='+',s=50,zorder=1)
                            counter[5] = counter[5] +1   
                        
        
        for i in range(len(counter)):
            counter[i] = int(counter[i] / 2)

        print('\n')
        for idx2,name in enumerate(['spitzer class I', 'spitzer class II', 'spitzer no class', 'match class I', 'match class II', 'match no class','almagal source']):
            print(f'{name}: {counter[idx2]}')
                    
        print('\n')
    else:
    
    # MARK: ALMA ELSE
        
        class_I = []
        class_II = []

        #spitzer 
        for i in range(2):
            if i == 1:
                class_I = []
                class_II = []

            for source in spitzer_sources:
                if source['p1_class_I'] == 'True':
                    axs[i].scatter(source['l']*u.deg-center[0],source['b']*u.deg-center[1],color='blue',marker='x',s=50,zorder=4)
                    class_I.append(((source['l']*u.deg-center[0]).value,(source['b']*u.deg-center[1]).value))
                elif source['p1_class_II'] == 'True':
                    axs[i].scatter(source['l']*u.deg-center[0],source['b']*u.deg-center[1],color='red',marker='x',s=50,zorder=3)
                    class_II.append(((source['l']*u.deg-center[0]).value,(source['b']*u.deg-center[1]).value))
                else:
                    axs[i].scatter(source['l']*u.deg-center[0],source['b']*u.deg-center[1],color='lightgrey',marker='x',s=50,zorder=2)
            

            axs[i].scatter(alma_sources['GLON']*u.deg-center[0],alma_sources['GLAT']*u.deg-center[1],color='grey',marker='+',s=50,zorder=1)
        
        print('\n')
    
        print(f'spitzer class I: {len(class_I)}')
        print(f'spitzer class II: {len(class_II)}')
        print(f'spitzer no class: {len(spitzer_sources) - len(class_I) - len(class_II)}')    

        print('\n')

    class_I_bol, class_II_bol, false_dist = False, False, False

    # MARK: Centroid

    if len(class_I) > 3:
        centroid_class_I = shapely.centroid(MultiPoint(class_I))
        class_I_bol = True
    if len(class_II) > 3:
        centroid_class_II = shapely.centroid(MultiPoint(class_II))
        class_II_bol = True
    
    if class_I_bol and class_II_bol:
        centroid_total = shapely.centroid(LineString([(centroid_class_I.x,centroid_class_I.y),(centroid_class_II.x,centroid_class_II.y)]))
    elif class_I_bol:
        centroid_total = centroid_class_I
    elif class_II_bol:
        centroid_total = centroid_class_II
    else:
        centroid_total = Point(0,0)
        false_dist = True
    
    if not false_dist:
        for i in range(2):
            axs[i].scatter(centroid_total.x,centroid_total.y,color='magenta',marker='D',s=30)
            axs[i].plot([centroid_total.x,0],[centroid_total.y,0],color='magenta',linestyle='--',alpha=0.5)
            
    dist = np.sqrt(centroid_total.x**2+centroid_total.y**2)

    print('\n')
    print(f'distance to centroid: {dist} [deg]')

    axs[0].indicate_inset_zoom(axs[1],edgecolor='red')
    fig.suptitle(f'Region {region}',y=0.72)

    os.makedirs('../output/figures/'+cat+'/regions',exist_ok=True)
    os.makedirs('../output/figures/'+cat+'/regions/match',exist_ok=True)
    os.makedirs('../output/figures/'+cat+'/regions/spitzer',exist_ok=True)
    if alma_proceeding:
        fig.savefig('../output/figures/'+cat+'/regions/match/'+str(int(idx)+1)+'_'+region+'_match.pdf',format='pdf',dpi=300,bbox_inches='tight')
        return [dist,false_dist,counter[3],counter[4]]
    else:
        fig.savefig('../output/figures/'+cat+'/regions/spitzer/'+str(int(idx)+1)+'_'+region+'_spitzer.pdf',format='pdf',dpi=300,bbox_inches='tight')
        return [dist,false_dist,len(class_I),len(class_II)]







    



