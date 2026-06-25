from astropy.table import Table
from YSO_statistics import YSO_plots_compare, region_analyze, region_plot
import numpy as np
from matplotlib import pyplot as plt 

cat = input("\ninput catalog name (leave blank for glimpse_s07) ").lower()

if cat == '':
    cat = 'glimpse_s07'

alma_proceeding = (input("Do you want to work with ALMA data? (y/n): ").strip().lower() == "y")

#alma_proceeding = True
#cat = 'glimpse_s07'

if alma_proceeding:
    spitzer_sources_alma = Table.read("../output/tables/"+cat+"/spitzer_sources_YSO_PLOTS_ALMA.tbl", format="ascii.ipac")
    alma_regions_alma = Table.read('../output/tables/'+cat+'/alma_regions_YSO_PLOTS_ALMA.tbl', format='ascii.ipac')
    alma_sources_alma = Table.read('../output/tables/'+cat+'/alma_sources_YSO_PLOTS_ALMA.tbl', format='ascii.ipac')

spitzer_sources = Table.read("../output/tables/"+cat+"/spitzer_sources_YSO_PLOTS.tbl", format="ascii.ipac")
alma_regions = Table.read('../output/tables/'+cat+'/alma_regions_YSO_PLOTS.tbl', format='ascii.ipac')
alma_sources = Table.read('../output/tables/'+cat+'/alma_sources_YSO_PLOTS.tbl', format='ascii.ipac')

region_info = Table.read('../output/tables/alma_regions_coordinates.tbl',format='ascii.ipac')

if input('do you want to analyze the regions (y/n) ').lower() == 'y':
    if alma_proceeding:
        alma_regions_alma, region_id_total, region_id_YSO = region_analyze(spitzer_sources_alma,alma_regions_alma,cat,alma_proceeding)
    else:
        alma_regions, region_id_total, region_id_YSO = region_analyze(spitzer_sources,alma_regions,cat,alma_proceeding)

    if input('do you want to plot them? (y/n) ').lower('y') == 'y':
        n_regions = input('how many regions? (max 1017) (caution: long time to compile) ')
        data_region_id = []
        data_dist = []
        data_false_dist = []
        data_class_I = []
        data_class_II = []


        for idx,region in enumerate(region_id_YSO[:n_regions]):

            if alma_proceeding:
                data = region_plot(spitzer_sources[spitzer_sources['internal_id_01'] == region],alma_sources[alma_sources['CLUMP'] == region],alma_sources_alma[alma_sources_alma['CLUMP'] == region],cat,region,region_info[region_info['INTERNAL_ID'] == region],alma_proceeding,idx)
            else:
                data = region_plot(spitzer_sources[spitzer_sources['internal_id_01'] == region],alma_sources[alma_sources['CLUMP'] == region],alma_sources[alma_sources['CLUMP'] == region],cat,region,region_info[region_info['INTERNAL_ID'] == region],alma_proceeding,idx)

            data_region_id.append(region)
            data_dist.append(data[0])
            data_false_dist.append(data[1])
            data_class_I.append(data[2])
            data_class_II.append(data[3])

        tbl = Table()

        tbl['region_id'] = data_region_id
        tbl['distance'] = data_dist
        tbl['false_dist'] = data_false_dist
        tbl['class_I'] = data_class_I
        tbl['class_II'] = data_class_II

        if alma_proceeding:

            tbl.write(
                '../output/tables/'+cat+'/distance_match.tbl',
                format='ascii.ipac',
                overwrite=True
            )
        else:
            tbl.write(
                '../output/tables/'+cat+'/distance_spitzer.tbl',
                format='ascii.ipac',
                overwrite=True
            )
if alma_proceeding:
    if input('do you want to plot histogram of sources, L/M, n(H2) plots? (y/n) ').lower() == 'y':
        YSO_plots_compare(spitzer_sources,alma_sources,spitzer_sources_alma,cat)

if input('do you want to plot centroid distance vs YSO? (y/n) ').lower() == 'y':
    if alma_proceeding:
        tbl = Table.read('../output/tables/'+cat+'/distance_match.tbl',format='ascii.ipac')
        tbl = tbl[tbl['false_dist'] == 'False']
        fig, axs = plt.subplots(2,2,figsize=(3,3),gridspec_kw={
            'width_ratios': [4, 1],
            'height_ratios': [1, 4]
        })

        axs[0][1].axis('off')
        axs[1][0].hist2d(tbl['distance']*3600,tbl['class_I']+tbl['class_II'],cmap='plasma',bins=[np.array(range(0,45,5)),range(0,5,1)])

        axs[0][0].hist(tbl['distance']*3600,bins=np.array(range(0,450,25))/10)
        axs[1][1].hist(tbl['class_I']+tbl['class_II'],bins=range(0,5,1),orientation='horizontal')


        axs[0][0].set_xlim([0,40])
        axs[0][0].set_ylim([0,12])
        axs[1][1].set_ylim([0,4])
        axs[1][1].set_xlim([0,50])

        axs[1][0].set_xlim([0,40])
        axs[1][0].set_ylim([0,4])

        counts, _ = np.histogram(tbl['class_I']+tbl['class_II'],bins=range(0,np.max(tbl['class_I']+tbl['class_II'])+1,1))

        print(counts)

        axs[1][0].set_xlabel('distance from centroid [arcsec]')
        axs[1][0].set_ylabel('YSO detected (matched data)')

        axs[0][0].tick_params(axis='x', bottom=False, labelbottom=False)
        axs[1][1].tick_params(axis='y', left=False, labelleft=False)
        axs[1][1].tick_params(axis='y', which='minor',left=False, labelleft=False)

        fig.savefig('../output/figures/'+cat+'/distance_YSO_ALMA.pdf',format='pdf',dpi=300)
    else:
        tbl = Table.read('../output/tables/'+cat+'/distance_spitzer.tbl',format='ascii.ipac')
        tbl = tbl[tbl['false_dist'] == 'False']
        fig, axs = plt.subplots(2,2,figsize=(3,3),gridspec_kw={
            'width_ratios': [4, 1],
            'height_ratios': [1, 4]
        })

        axs[0][1].axis('off')
        axs[1][0].hist2d(tbl['distance']*3600,tbl['class_I']+tbl['class_II'],cmap='plasma',bins=[range(4,42,4),range(4,15,1)])

        axs[0][0].hist(tbl['distance']*3600,bins=range(4,42,2))
        axs[1][1].hist(tbl['class_I']+tbl['class_II'],orientation='horizontal',bins=range(4,15,1))


        axs[0][0].set_xlim([4,40])
        axs[0][0].set_ylim([0,14])
        

        axs[1][1].set_xlim([4,22])
        axs[1][1].set_ylim([4,14])

        axs[1][0].set_ylim([4,14])
        axs[1][0].set_xlim([4,40])

        counts, _ = np.histogram(tbl['class_I']+tbl['class_II'],bins=range(4,np.max(tbl['class_I']+tbl['class_II']),1))

        print(counts)
        
        axs[1][0].set_xlabel('distance from centroid [arcsec]')
        axs[1][0].set_ylabel('YSO detected (Spitzer data)')

        axs[0][0].tick_params(axis='x', bottom=False, labelbottom=False)
        axs[1][1].tick_params(axis='y', left=False, labelleft=False)
        axs[1][1].tick_params(axis='y', which='minor',left=False, labelleft=False)

        fig.savefig('../output/figures/'+cat+'/distance_YSO_spitzer.pdf',format='pdf',dpi=300)
        


