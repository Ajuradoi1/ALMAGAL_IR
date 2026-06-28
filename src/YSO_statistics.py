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


#MARK: Region analyzer

def region_analyze(spitzer_sources,spitzer_sources_alma,alma_regions,alma_sources,alma_sources_alma,cat):
    
    cols = ['pah_galaxies', 'agn_sources', 'sed_sources', 'pah_contaminated','p1_class_I', 'p1_class_II', 'p1_noclass', 'p2_class_I', 'p2_class_II', 'p2_noclass'] 

    count_60_class_I = []
    count_60_class_II = []
    count_20_class_0 = []
    count_20_class_I  = []
    count_20_class_II = []
    count_20_radio = []

    ratio_60_class_I = []
    ratio_60_class_II = []
    ratio_20_class_0 = []
    ratio_20_class_I  = []
    ratio_20_class_II = []
    ratio_20_radio = []

    count_debug = 0
    for col in cols:
        alma_regions[col] = 0
    
    for idx, region in enumerate(alma_regions):
        region_id = region['INTERNAL_ID']

        count_60_class_I.append(0)
        count_60_class_II.append(0)
        count_20_class_0.append(0)
        count_20_class_I.append(0)
        count_20_class_II.append(0)
        count_20_radio.append(0)


        count_60_class_I[-1] = len(spitzer_sources[(spitzer_sources['p1_class_I'] == 'True') & (spitzer_sources['internal_id_01'] == region_id)])
        count_60_class_II[-1] =  len(spitzer_sources[(spitzer_sources['p1_class_II'] == 'True') & (spitzer_sources['internal_id_01'] == region_id)])

        

        count_20_class_I[-1] = len(spitzer_sources_alma[(spitzer_sources_alma['p1_class_I'] == 'True') & (spitzer_sources_alma['internal_id_01'] == region_id)]) + len(spitzer_sources_alma[(spitzer_sources_alma['p2_class_I'] == 'True') & (spitzer_sources_alma['internal_id_01'] == region_id)])
        count_20_class_II[-1] =  len(spitzer_sources_alma[(spitzer_sources_alma['p1_class_II'] == 'True') & (spitzer_sources_alma['internal_id_01'] == region_id)]) + len(spitzer_sources_alma[(spitzer_sources_alma['p2_class_II'] == 'True') & (spitzer_sources_alma['internal_id_01'] == region_id)])
            
        count_20_class_0[-1] = len(alma_sources[(alma_sources['CLUMP'] == region_id) & (alma_sources['HAS_SPITZER'] == 'False') & (alma_sources['RADIO_MATCH'] == 0)]) 
        count_20_radio[-1] = len(alma_sources[(alma_sources['CLUMP'] == region_id) & (alma_sources['HAS_SPITZER'] == 'False') & (alma_sources['RADIO_MATCH'] == 1)])
        # Compute % of Class I and Class II for 60'' (Total spitzer data)

        aux_val = (count_60_class_I[-1]+count_60_class_II[-1])

        

        if aux_val != 0:

            ratio_60_class_I.append(count_60_class_I[-1] / aux_val)
            ratio_60_class_II.append(count_60_class_II[-1] / aux_val)
        else:
            ratio_60_class_I.append(0)
            ratio_60_class_II.append(0)

        # Compute % of Class I and Class II for 20'' (ALMAGAL)

        aux_val2 =  (count_20_class_0[-1]+count_20_class_I[-1]+count_20_class_II[-1]+count_20_radio[-1])

        if aux_val2 != 0:
            ratio_20_class_0.append((count_20_class_0[-1]) / aux_val2)
            ratio_20_class_I.append((count_20_class_I[-1]) /aux_val2)
            ratio_20_class_II.append((count_20_class_II[-1]) / aux_val2)
            ratio_20_radio.append(count_20_radio[-1] / aux_val2)
        else:
            ratio_20_class_0.append(0)
            ratio_20_class_I.append(0)
            ratio_20_class_II.append(0)
            ratio_20_radio.append(0)
            

    alma_regions['class_I_60'] = count_60_class_I
    alma_regions['class_II_60'] = count_60_class_II

    alma_regions['class_0_20'] = count_20_class_0
    alma_regions['class_I_20'] = count_20_class_I
    alma_regions['class_II_20'] = count_20_class_II
    alma_regions['radio_20'] = count_20_radio

    alma_regions['rat_class_I_60'] = ratio_60_class_I
    alma_regions['rat_class_II_60'] = ratio_60_class_II

    alma_regions['rat_class_0_20'] = ratio_20_class_0
    alma_regions['rat_class_I_20'] = ratio_20_class_I
    alma_regions['rat_class_II_20'] = ratio_20_class_II
    alma_regions['rat_radio_20'] = ratio_20_radio


    alma_regions.write('../output/tables/'+cat+'/alma_regions_YSO.tbl',format='ascii.ipac',overwrite=True)

    count_dict = {}

    '''for idx,region in enumerate(alma_regions):
        
        count_60 = 0
        count_20 = 0
        for col in [['class_I_60','class_II_60']],:
            count_60 = count_60 + region[col[0]]+region[col[1]]
        
        for col in [['class_I_20','class_II_20']]:
                count_20 = count_20 + region[col[0]]+region[col[1]]'''
    
    count_60 = np.array(alma_regions['class_I_60'])+np.array(alma_regions['class_I_60'])
    count_20 = np.array(alma_regions['class_I_20'])+np.array(alma_regions['class_II_20'])


    '''total_tol = sorted(count_dict.items(),key=lambda item: item[1][1],reverse=True)
    total_YSO = sorted(count_dict.items(),key=lambda item: item[1][2],reverse=True)

    total_id = []
    yso_id = []

    for element in total_tol:
        total_id.append(element[1][0])
    for element in total_YSO:
        yso_id.append(element[1][0])

    '''

    tbl = Table()

    tbl['region_id'] = alma_regions['INTERNAL_ID']
    tbl['class_60'] = count_60
    tbl['class_20'] = count_20
    
    #print(f'Region with most sources: {total_id[0]}\n')
    #print(f'Region with most YSO matched: {yso_id[0]}\n')

    tbl1 = tbl.copy()

    tbl1.sort('class_60')
    tbl1.reverse()

    tbl2 = tbl.copy()

    tbl2.sort('class_20')
    tbl2.reverse()

    return alma_regions, np.array(tbl1['region_id']), np.array(tbl2['region_id'])


def YSO_plots_compare(spitzer_sources,alma_regions,spitzer_sources_alma,cat):


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

    fig, axs = plt.subplots(2,1,figsize=(6,4),tight_layout=True,gridspec_kw={'height_ratios': [7, 3]},sharex=True)

    axs[0].tick_params(axis='x', bottom=False,labelbottom=False)

    

    fig.subplots_adjust(bottom=0.25)
    fig.subplots_adjust(right=0.8)
    #fig.subplots_adjust(left=0.7999)
    
    for idx, table in enumerate([spitzer_sources,spitzer_sources_alma]):

        counts = []

        for col in cols:
            
            counts.append(len(table[col][table[col] == 'True']))


        x = np.arange(len(counts))

        if idx == 0:
            bars = axs[0].bar(x, counts,color='lightgrey')
            axs[0].axhline(sum(counts),color='lightgrey',linestyle='--',alpha=0.5)
            axs[0].text(7,sum(counts),f"{sum(counts)} sources",ha='left',va='bottom',fontsize=9)
        else:
            
            bars = axs[0].bar(x, counts,hatch='///',facecolor='none',edgecolor='magenta')
            axs[0].axhline(sum(counts),color='magenta',linestyle='--',alpha=0.5)
            axs[0].text(7,sum(counts),f"{sum(counts)} sources",ha='left',va='bottom',fontsize=9)

        

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
            axs[0].set_ylim(1,10 ** len(str(sum(counts))))

    axs[1].set_xticks(x)
    axs[1].set_xticklabels(['PAH Gal.','AGN','Shock', 'PAH Contam.', 'Class I', 'Class II', 'No Class'], rotation=45, ha="right", fontsize=10)

    

    counts = []

    cols_2 = ["p2_class_I", "p2_class_II", "p2_noclass"]
    for col in cols_2:
        counts.append(len(table[col][table[col] == 'True']))
    
    x = np.arange(len(counts)) + 4

    bars = axs[1].bar(x, counts,hatch='///',facecolor='none',edgecolor='magenta')
    axs[1].axhline(sum(counts),color='magenta',linestyle='--',alpha=0.5)
    axs[1].text(7,sum(counts),f"{sum(counts)} sources",ha='left',va='bottom',fontsize=9)

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

    axs[0].set_yscale('log')
    axs[1].set_yscale('log')
    axs[1].set_ylim(1,10 ** len(str(sum(counts))))

    axs[0].set_ylabel('Phase 1 ')
    axs[1].set_ylabel('Phase 2 ')

    fig.savefig("../output/figures/"+cat+"/YSO_classification_histogram_"+cat+".pdf",format='pdf',dpi=300)
    plt.close()
    
    # MARK: L/M, L, M, T, DIST, N(H_2)

    class_I_60 = []
    class_II_60 = []
    class_0_20 = []
    class_I_20 = []
    class_II_20 = []
    radio_20 = []

    SURFD_to_nH2 = (
        
        6.022e23  # Avogadro's constant (N_A) [particles/mol]
        / 2.016  # molecular mass of n(H2) [g/mol]
    )

    for region in alma_regions:

        if (
            region['MASS'] == -999 or
            region['LUM'] == -999 or
            region['TBOL'] == -999 or 
            region['DIST_AG'] == -999 or
            region['NODEC_SURFD'] == -999
        ):
            continue
        
        L_M = np.log10(region['LUM'] / region['MASS'])
        L = np.log10(region['LUM'])
        M = np.log10(region['MASS'])
        dist = np.log10(region['DIST_AG'])
        T = np.log10(region['TBOL'])
        NH_2 = np.log10(region['NODEC_SURFD']*SURFD_to_nH2)


        class_I_60.append([L_M,L,M,dist,T,NH_2, region['rat_class_I_60']*100])
        class_II_60.append([L_M,L,M,dist,T,NH_2, region['rat_class_II_60']*100])
        class_0_20.append([L_M,L,M,dist,T,NH_2, region['rat_class_0_20']*100])
        class_I_20.append([L_M,L,M,dist,T,NH_2, region['rat_class_I_20']*100])
        class_II_20.append([L_M,L,M,dist,T,NH_2, region['rat_class_II_20']*100])
        radio_20.append([L_M,L,M,dist,T,NH_2, region['rat_radio_20']*100])

    xlabels = [r'log L/M [L$_\odot$/M$_\odot$]', r'log L [L$_\odot$]', r'log M [M$_\odot$]', 'log distance [pc]', r'log T$_{\text{bol}}$ [K]', r'log N(H$_2$) [cm$^{-2}$]']

    for idx, plot in enumerate(['L_M','L','M','dist','T','NH_2']):

        fig, axs = plt.subplots(2,1,figsize=(5,4),tight_layout=True,sharex=True)

        for idx2,data_set in enumerate([[class_I_60,class_II_60,class_0_20, radio_20],[class_I_20,class_II_20,class_0_20, radio_20]]):

            data_1 = [[x[idx] for x in data_set[0]],[x[-1] for x in data_set[0]]]
            data_2 = [[x[idx] for x in data_set[1]],[x[-1] for x in data_set[1]]]
            data_3 = [[x[idx] for x in data_set[2]],[x[-1] for x in data_set[2]]]
            data_4 = [[x[idx] for x in data_set[3]],[x[-1] for x in data_set[3]]]

            x1, y1 = np.array(data_1[0]),np.array(data_1[1])
            x2, y2 = np.array(data_2[0]), np.array(data_2[1])
            x3, y3 = np.array(data_3[0]), np.array(data_3[1])
            x4, y4 = np.array(data_4[0]), np.array(data_4[1])
            

            if idx2 == 0:
                min = round(10*np.min([np.min(x1),np.min(x2),np.min(x3),np.min(x4)]))/10
                max = round(10*np.max([np.max(x1),np.max(x2),np.max(x3),np.max(x4)]))/10
                n_bins = 12
                step = (max-min)/(n_bins)
                bins = np.arange(min, max+step, step)

            

            groups1, groups2, groups3,groups4 = [], [], [], []
            pos1, pos2, pos3, pos4 = [], [], [], []

            #offset = step/50  # spacing between boxplots

            offset = 0

            for i in range(len(bins) - 1):
                mask1 = (x1 >= bins[i]) & (x1 < bins[i+1])
                mask2 = (x2 >= bins[i]) & (x2 < bins[i+1])
                mask3 = (x3 >= bins[i]) & (x3 < bins[i+1])
                mask4 = (x4 >= bins[i]) & (x4 < bins[i+1])

                center = (bins[i] + bins[i+1]) / 2

                if np.sum(mask1) > 0:
                    pos1.append(center)
                    groups1.append(y1[mask1])
                    
                if np.sum(mask2) > 0:
                    pos2.append(center + offset)
                    groups2.append(y2[mask2])
                if np.sum(mask3) > 0:
                    pos3.append(center - offset)       
                    groups3.append(y3[mask3])
                if np.sum(mask4) > 0:
                    pos4.append(center - 2*offset)       
                    groups4.append(y4[mask4])


            # Class I

            bp1 = axs[idx2].boxplot(
                groups1,
                positions=pos1,
                widths=step,
                patch_artist=True,
                showfliers=False,
                zorder=2
            )
            axs[idx2].scatter(x1,y1,color='blue',s=10,alpha=0.6,zorder=20)
            for box in bp1['boxes']:
                box.set(facecolor='blue', alpha=0.5, edgecolor='black',linewidth=2,zorder=0)

            # Class II

            bp2 = axs[idx2].boxplot(
                groups2,
                positions=pos2,
                widths=step,
                patch_artist=True,
                showfliers=False,
                zorder=3
            )
            axs[idx2].scatter(x2,y2,color='red',s=10,alpha=0.6,zorder=20)
            for box in bp2['boxes']:
                box.set(facecolor='red', alpha=0.4, edgecolor='red',hatch='///',linewidth=1,zorder=1)


        # Class 0

        bp0 = axs[1].boxplot(
            groups3,
            positions=pos3,
            widths=step,
            patch_artist=True,
            showfliers=False,
            zorder=1
        )
        
        for box in bp0['boxes']:
            box.set(facecolor='limegreen', alpha=0.5,edgecolor='black')
        axs[1].scatter(x3,y3,color='darkgreen',s=10,alpha=0.3,zorder=20)

        bp3 = axs[1].boxplot(
            groups4,
            positions=pos4,
            widths=step,
            patch_artist=True,
            showfliers=False,
            zorder=0
        )
        axs[1].scatter(x4,y4,color='orange',s=10,alpha=0.6,zorder=20)
        for box in bp3['boxes']:
            box.set(facecolor='gold', alpha=0.4,edgecolor='black',hatch='///')


        axs[0].set_xlim(min,max)
        axs[0].set_ylim([-5,105])
        axs[1].set_ylim([-5,105])

        axs[1].set_xlabel(f'{xlabels[idx]}')
        axs[0].set_ylabel(r'% of YSO [Spitzer]')
        axs[1].set_ylabel(r'% of YSO [ALMA]')

        ticks = bins[::2]

        axs[1].set_xticks(ticks, [f"{t:.2f}" for t in ticks])
        if idx2 == 1:
            fig.savefig('../output/figures/'+cat+'/results_'+plot+'.pdf',format='pdf',dpi=300)
            plt.show()
            plt.close()


#MARK: region_plot

def region_plot(spitzer_sources,alma_sources,alma_sources_alma,cat,region,region_info,mode,idx): # to fix

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

    alma_proceeding = True

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
                        if (source['HAS_SPITZER'] =='False') and (source['RADIO_MATCH'] == 0):
                            axs[i].scatter(table['GLON']*u.deg-center[0],table['GLAT']*u.deg-center[1],color='limegreen',marker='+',s=50,zorder=2)
                            counter[6] = counter[6] +1
                        elif (source['HAS_SPITZER'] == 'False') and (source['RADIO_MATCH'] == 1):
                            axs[i].scatter(table['GLON']*u.deg-center[0],table['GLAT']*u.deg-center[1],color='gold',marker='+',s=50,zorder=3)
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
    
    dist_class_I,dist_class_II = [],[]
    if not false_dist:
        for i in range(2):
            axs[i].scatter(centroid_total.x,centroid_total.y,color='magenta',marker='D',s=30)
            axs[i].plot([centroid_total.x,0],[centroid_total.y,0],color='magenta',linestyle='--',alpha=0.5)
        
        for source in spitzer_sources:
            pit_dist = np.sqrt((centroid_total.x-(source['l']*u.deg-center[0]).value)**2+(centroid_total.y-(source['b']*u.deg-center[1]).value)**2)

            if source['p1_class_I'] == 'True':
                dist_class_I.append(pit_dist)
            elif source['p1_class_II'] == 'True':
                dist_class_II.append(pit_dist)


        if len(dist_class_I) != 0:
            dist_class_I_np = np.mean(np.array(dist_class_I))

        if len(dist_class_II) != 0:
            dist_class_II_np = np.mean(np.array(dist_class_II))

        if len(dist_class_I) == 0:
            dist_class_I = -999
        else:
            dist_class_I = dist_class_I_np

        if len(dist_class_II) == 0:
            dist_class_II = -999
        else:
            dist_class_II = dist_class_II_np 
          
           


    else:
        dist_class_I,dist_class_II  = [-999,-999]


    dist_class_0,dist_radio = [],[]
    dist_class_I_alma,dist_class_II_alma = [],[]

    counter_class_0 = 0
    counter_radio = 0
    counter_class_I_alma = 0
    counter_class_II_alma = 0


    if not false_dist:
        for i,sources in enumerate([alma_sources,alma_sources_alma]):
            for source in sources:
                pit_dist = np.sqrt((centroid_total.x-(source['GLON']*u.deg-center[0]).value)**2+(centroid_total.y-(source['GLAT']*u.deg-center[1]).value)**2)
                
                if i == 0:
                    
                    if (source['HAS_SPITZER'] == 'False') and (source['RADIO_MATCH'] == 0):
                        dist_class_0.append(pit_dist)
                        counter_class_0 = counter_class_0 +1
                    elif (source['HAS_SPITZER'] == 'False') and (source['RADIO_MATCH'] == 1):
                        dist_radio.append(pit_dist)
                        counter_radio = counter_radio +1

                else:
                    if ((source['p1_class_I'] == 'True') or (source['p2_class_I'] == 'True')):
                        
                        dist_class_I_alma.append(pit_dist)
                        counter_class_I_alma = counter_class_I_alma + 1

                    elif (source['p1_class_II'] == 'True') or (source['p2_class_II'] == 'True'):
                        
                        dist_class_II_alma.append(pit_dist)
                        counter_class_II_alma = counter_class_II_alma + 1 



        if len(dist_class_0) != 0:
            dist_class_0_np = np.mean(np.array(dist_class_0))

        if len(dist_radio) != 0:
            dist_radio_np = np.mean(np.array(dist_radio))

        if len(dist_class_I_alma) != 0:
            dist_class_I_alma_np = np.mean(np.array(dist_class_I_alma))

        if len(dist_class_II_alma) != 0:
            dist_class_II_np = np.mean(np.array(dist_class_II_alma))

        if len(dist_class_0) == 0:
            dist_class_0 = -999
        else:
            dist_class_0 = dist_class_0_np

        if len(dist_radio) == 0:
            dist_radio = -999
        else:
            dist_radio = dist_radio_np

        if len(dist_class_I_alma) == 0:
            dist_class_I_alma = -999
        else:
            dist_class_I_alma = dist_class_I_alma_np 
        
        if len(dist_class_II_alma) == 0:
            dist_class_II_alma = -999
        else:
            dist_class_II_alma = dist_class_II_np 
          
    else:
        dist_class_0,dist_radio  = [-999,-999]    
        dist_class_I_alma,dist_class_II_alma = [-999,-999]
    


    axs[0].indicate_inset_zoom(axs[1],edgecolor='red')
    fig.suptitle(f'Region {region}',y=0.72)

    os.makedirs('../output/figures/'+cat+'/regions',exist_ok=True)
    os.makedirs('../output/figures/'+cat+'/regions/match',exist_ok=True)
    os.makedirs('../output/figures/'+cat+'/regions/spitzer',exist_ok=True)
    if mode == 'match':
        fig.savefig('../output/figures/'+cat+'/regions/match/'+str(int(idx)+1)+'_'+region+'_match.pdf',format='pdf',dpi=300,bbox_inches='tight')
        plt.close()
    else:
        fig.savefig('../output/figures/'+cat+'/regions/spitzer/'+str(int(idx)+1)+'_'+region+'_spitzer.pdf',format='pdf',dpi=300,bbox_inches='tight')
        plt.close()

    return [[dist_class_I,dist_class_II,dist_class_0,dist_radio,dist_class_I_alma,dist_class_II_alma],[len(class_I),len(class_II),counter_class_0,counter_radio,counter_class_I_alma,counter_class_II_alma],false_dist]


# MARK: distribution_plots

def distribution_plots(n_regions,mode,cat):

    tbl = Table.read('../output/tables/'+cat+'/'+str(n_regions)+'_distance_'+str(mode)+'.tbl',format='ascii.ipac')

    tbl = tbl[tbl['false_dist'] == 'False']

    

    aux = np.array(tbl['class_I'])+np.array(tbl['class_II'])

    aux[aux == 0] = 1

    count_I = np.array(tbl['class_I']) / aux  
    count_II = np.array(tbl['class_II']) / aux

    aux2 = np.array(tbl['class_0'])+np.array(tbl['radio'])+np.array(tbl['class_I_alma'])+np.array(tbl['class_II_alma'])

    aux2[aux2 == 0] = 1  
    
    count_0 = np.array(tbl['class_0']) / aux2  
    count_radio = np.array(tbl['radio']) / aux2
    count_I_alma = np.array(tbl['class_I_alma']) / aux2  
    count_II_alma = np.array(tbl['class_II_alma']) / aux2


    for idx, count in enumerate(count_I):
        if np.isnan(count):
            count_I[idx] = 0

    for idx, count in enumerate(count_II):
        if np.isnan(count):
            count_II[idx] = 0

    for idx, count in enumerate(count_0):
        if np.isnan(count):
            count_0[idx] = 0

    for idx, count in enumerate(count_radio):
        if np.isnan(count):
            count_radio[idx] = 0

    for idx, count in enumerate(count_I_alma):
        if np.isnan(count):
            count_I_alma[idx] = 0

    for idx, count in enumerate(count_II_alma):
        if np.isnan(count):
            count_II_alma[idx] = 0



    tbl['rat_class_I'] = count_I
    tbl['rat_class_II'] = count_II

    tbl['rat_class_0'] = count_0
    tbl['rat_radio'] = count_radio
    tbl['rat_class_I_alma'] = count_I_alma
    tbl['rat_class_II_alma'] = count_II_alma

    tbl.write('../output/tables/glimpse_s07/final_plot.tbl',format='ascii.ipac',overwrite=True)

    fig, axs = plt.subplots(2,1,figsize=(5,4),tight_layout=True)

    color_list = [['blue','red'], ['lightgrey','darkgrey']]

    for idx, data_set in enumerate([
        [['dist_class_I','dist_class_II','dist_class_I_alma','dist_class_II_alma'],['rat_class_I','rat_class_II','rat_class_I_alma','rat_class_II_alma']],
        [['dist_class_I_alma','dist_class_II_alma','dist_class_0','dist_class_radio'],['rat_class_I_alma','rat_class_II_alma','rat_class_0','rat_radio',]]
        ]):

        if idx == 0:
            data = [
                [tbl[data_set[0][0]][tbl[data_set[0][0]] != -999], tbl[data_set[1][0]][tbl[data_set[0][0]] != -999]],
                [tbl[data_set[0][1]][tbl[data_set[0][1]] != -999], tbl[data_set[1][1]][tbl[data_set[0][1]] != -999]],
                [tbl[data_set[0][2]][tbl[data_set[0][2]] != -999], tbl[data_set[1][2]][tbl[data_set[0][2]] != -999]],
                [tbl[data_set[0][3]][tbl[data_set[0][3]] != -999], tbl[data_set[1][3]][tbl[data_set[0][3]] != -999]]
            ]
            
            x1, y1 = np.array(data[0][0])*3600,np.array(data[0][1])*100
            x2, y2 = np.array(data[1][0])*3600,np.array(data[1][1]) *100
            x3, y3 = np.array(data[2][0])*3600,np.array(data[2][1])*100
            x4, y4 = np.array(data[3][0])*3600, np.array(data[3][1]) *100
        else:
            data = [
                [tbl[data_set[0][0]][tbl[data_set[0][0]] != -999], tbl[data_set[1][0]][tbl[data_set[0][0]] != -999]],
                [tbl[data_set[0][1]][tbl[data_set[0][1]] != -999], tbl[data_set[1][1]][tbl[data_set[0][1]] != -999]],
                [tbl[data_set[0][2]][tbl[data_set[0][2]] != -999], tbl[data_set[1][2]][tbl[data_set[0][2]] != -999]],
                [tbl[data_set[0][3]][tbl[data_set[0][3]] != -999], tbl[data_set[1][3]][tbl[data_set[0][3]] != -999]]
            ]

            x1, y1 = np.array(data[0][0])*3600,np.array(data[0][1])*100
            x2, y2 = np.array(data[1][0])*3600,np.array(data[1][1]) *100
            x3, y3 = np.array(data[2][0])*3600,np.array(data[2][1])*100
            x4, y4 = np.array(data[3][0])*3600, np.array(data[3][1]) *100

        if idx == 0:
            min = np.min([np.min(x1),np.min(x2),np.min(x3),np.min(x4)])
            max = np.max([np.max(x1),np.max(x2),np.max(x3),np.max(x4)])
        if idx == 1:
            min_2 = np.min([np.min(x1),np.min(x2),np.min(x3),np.min(x4)])
            max_2 = np.max([np.max(x1),np.max(x2),np.max(x3),np.max(x4)])

            min = np.min([min,min_2])
            max = np.min([max,max_2])
        print(f'{min}: {max}')
        n_bins = 15
        step = (max-min)/(n_bins+1)
        bins = np.arange(min, max+step, step)

        groups1, groups2,groups3, groups4 = [], [], [], []
        pos1, pos2, pos3, pos4 = [], [], [], []

        # offset = step/50  # spacing between boxplots

        for i in range(len(bins) - 1):
            mask1 = (x1 >= bins[i]) & (x1 < bins[i+1])
            mask2 = (x2 >= bins[i]) & (x2 < bins[i+1])
            mask3 = (x3 >= bins[i]) & (x3 < bins[i+1])
            mask4 = (x4 >= bins[i]) & (x4 < bins[i+1])

            center = (bins[i] + bins[i+1]) / 2

            if np.sum(mask1) > 0:
                pos1.append(center)
                groups1.append(y1[mask1])
            if np.sum(mask2) > 0:
                pos2.append(center)
                groups2.append(y2[mask2])
            if np.sum(mask3) > 0:
                pos3.append(center)
                groups3.append(y3[mask3])
            if np.sum(mask4) > 0:
                pos4.append(center)
                groups4.append(y4[mask4])

        
        
        if idx == 0:

            # Class I 

            bp1 = axs[0].boxplot(
                groups1,
                positions=pos1,
                widths=step,
                patch_artist=True,
                showfliers=False,
                zorder=2
            )

            for box in bp1['boxes']:
                box.set(facecolor='blue', alpha=0.5, edgecolor='black',linewidth=2,zorder=0)

            axs[0].scatter(x1,y1,color='blue',s=10,alpha=0.5)
            # Class II
            
            bp2 = axs[0].boxplot(
                groups2,
                positions=pos2,
                widths=step,
                patch_artist=True,
                showfliers=False,
                zorder=3
            )
            
            for box in bp2['boxes']:
                box.set(facecolor='red', alpha=0.5, edgecolor=color_list[idx][0],hatch='///',linewidth=1)
            axs[0].scatter(x2,y2,color='red',s=10,alpha=0.5)
        if idx == 1:
            
            # Class I ALMA

            bp1 = axs[1].boxplot(
                groups1,
                positions=pos1,
                widths=step,
                patch_artist=True,
                showfliers=False,
                zorder=2
            )

            for box in bp1['boxes']:
                box.set(facecolor='blue', alpha=0.5, edgecolor='black',linewidth=2,zorder=0)
            axs[1].scatter(x1,y1,color='blue',s=10,alpha=0.5)
            
            # Class II ALMA

            bp2 = axs[1].boxplot(
                groups2,
                positions=pos2,
                widths=step,
                patch_artist=True,
                showfliers=False,
                zorder=3
            )
            
            for box in bp2['boxes']:
                box.set(facecolor='red', alpha=0.5, edgecolor=color_list[idx][0],hatch='///',linewidth=1)
            axs[1].scatter(x2,y2,color='red',s=10,alpha=0.5)
            
            # Class 0

            bp3 = axs[1].boxplot(
                groups3,
                positions=pos3,
                widths=step,
                patch_artist=True,
                showfliers=False,
                zorder=2
            )
            axs[1].scatter(x3,y3,color='limegreen',s=10,alpha=0.5)
            for box in bp3['boxes']:
                box.set(facecolor='limegreen', alpha=0.5, edgecolor='black',linewidth=2,zorder=0)

            # radio

            bp4 = axs[1].boxplot(
                groups4,
                positions=pos4,
                widths=step,
                patch_artist=True,
                showfliers=False,
                zorder=3
            )
            axs[1].scatter(x4,y4,color='gold',s=10,alpha=0.5)
            for box in bp4['boxes']:
                box.set(facecolor='gold', alpha=0.5, edgecolor='black',hatch='///',linewidth=1,zorder=2)

        axs[idx].set_xlim(min,max)
        axs[0].set_ylim([-5,105])
        axs[1].set_ylim([-5,105])
        
        axs[1].set_xlabel('mean distance to centroid [arcsec]')

        if idx == 0:
            axs[0].set_ylabel(r'% of YSO [Spitzer]')
        else:
            axs[1].set_ylabel(r'% sources [ALMA]')
        
        _, connectors = axs[0].indicate_inset_zoom(axs[1],edgecolor='k',linestyle='--',linewidth=2)

        

        for c in connectors:
            c.set_visible(False)

        ticks = bins[::3]

        axs[idx].set_xticks(ticks, [f"{t:.2f}" for t in ticks])

    fig.savefig('../output/figures/'+cat+'/results_distance_YSO.pdf',format='pdf',dpi=300)


    



