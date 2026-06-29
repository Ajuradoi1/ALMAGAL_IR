# ALMAGAL_IR_PLOT. Need to run ALMAGAL_IR_YSO for tables needed to run. 2 runs needed. 1 ALMA data and other 1 SPITZER
# Author: Álex Jurado i López
# github repo: https://github.com/Ajuradoi1/ALMAGAL_IR


from astropy.table import Table
from YSO_statistics import YSO_plots_compare, region_analyze, region_plot, distribution_plots
import numpy as np
from matplotlib import pyplot as plt 

'''cat = input("\ninput catalog name (leave blank for glimpse_s07) ").lower()

if cat == '':
    cat = 'glimpse_s07'''

#alma_proceeding = (input("Do you want to work with ALMA data? (y/n): ").strip().lower() == "y")

cat = 'glimpse_s07'

spitzer_sources = Table.read("../output/tables/"+cat+"/spitzer_sources_to_analyze.tbl", format="ascii.ipac")
alma_regions = Table.read('../output/tables/'+cat+'/alma_regions_to_analyze.tbl', format='ascii.ipac')
alma_sources = Table.read('../output/tables/'+cat+'/alma_sources_to_analyze.tbl', format='ascii.ipac')

spitzer_sources_alma = Table.read("../output/tables/"+cat+"/spitzer_sources_to_analyze_ALMA.tbl", format="ascii.ipac")
alma_regions_alma = Table.read('../output/tables/'+cat+'/alma_regions_to_analyze_ALMA.tbl', format='ascii.ipac')
alma_sources_alma = Table.read('../output/tables/'+cat+'/alma_sources_to_analyze_ALMA.tbl', format='ascii.ipac')

alma_proceeding = True

region_info = Table.read('../output/tables/alma_regions_coordinates.tbl',format='ascii.ipac')

print('\n\nanalyzing regions...\n\n')

alma_regions, region_id_total, region_id_YSO = region_analyze(spitzer_sources,spitzer_sources_alma,alma_regions,alma_sources,alma_sources_alma,cat)

print(region_id_total[:3])
print('\n\n')
print(region_id_YSO[:3])

if input('do you want to plot histogram of sources, L/M, n(H2) plots? (y/n) ').lower() == 'y':

    print('\n\ngenerating plots...\n\n')
    YSO_plots_compare(spitzer_sources,alma_regions,spitzer_sources_alma,cat)

skip = False

if input('do you want to plot the regions distribution? (y/n) ').lower() == 'y':

    mode = '0'
    while (mode != 'match' and mode!= 'total'): 
        mode = input('choose sort (match/total): ').lower()

    n_regions = int(input('how many regions? (max 1017) (caution: long time to compile) '))
    
    data_region_id = []

    data_dist_class_I = []
    data_dist_class_II = []
    data_dist_class_0 = []
    data_dist_radio = []
    data_dist_class_I_alma = []
    data_dist_class_II_alma = []
    
    data_class_I = []
    data_class_II = []
    data_class_0 = []
    data_radio = []
    data_class_I_alma = []
    data_class_II_alma = []
    
    data_false_dist = []

    if mode == 'total':
        region_ids = region_id_total
    else:
        region_ids = region_id_YSO

    for idx,region in enumerate(region_ids[:n_regions]):

        data = region_plot(spitzer_sources[spitzer_sources['internal_id_01'] == region],alma_sources[alma_sources['CLUMP'] == region],alma_sources_alma[alma_sources_alma['CLUMP'] == region],cat,region,region_info[region_info['INTERNAL_ID'] == region],mode,idx)

        data_region_id.append(region)
        
        data_dist_class_I.append(data[0][0])
        data_dist_class_II.append(data[0][1])
        data_dist_class_0.append(data[0][2])
        data_dist_radio.append(data[0][3])
        data_dist_class_I_alma.append(data[0][4])
        data_dist_class_II_alma.append(data[0][5])
        
        data_class_I.append(data[1][0])
        data_class_II.append(data[1][1])
        data_class_0.append(data[1][2])
        data_radio.append(data[1][3])
        data_class_I_alma.append(data[1][4])
        data_class_II_alma.append(data[1][5])

        data_false_dist.append(data[-1])

        if idx == round(len(region_ids[:n_regions])/4):
            print('\n\n25"%" completed!\n\n')

        if idx == round(len(region_ids[:n_regions])/2):
            print('\n\n50"%" completed!\n\n')

        if idx == round(3*len(region_ids[:n_regions])/4):
            print('\n\n75"%" completed!\n\n')

    print('\n\n100"%" completed!\n\n')

    tbl = Table()

    tbl['region_id'] = data_region_id
    
    tbl['dist_class_I'] = data_dist_class_I
    tbl['dist_class_II'] = data_dist_class_II
    tbl['dist_class_0'] = data_dist_class_0
    tbl['dist_class_radio'] = data_dist_radio
    tbl['dist_class_I_alma'] = data_dist_class_I_alma
    tbl['dist_class_II_alma'] = data_dist_class_II_alma
    
    tbl['class_I'] = data_class_I
    tbl['class_II'] = data_class_II
    tbl['class_0'] = data_class_0
    tbl['radio'] = data_radio
    tbl['class_I_alma'] = data_class_I_alma
    tbl['class_II_alma'] = data_class_II_alma

    tbl['false_dist'] = data_false_dist
    # [dist_class_I,dist_class_II,counter[3],counter[4],false_dist]

    print(tbl)

    if mode == 'match':

        tbl.write(
            '../output/tables/'+cat+'/'+str(n_regions)+'_distance_match.tbl',
            format='ascii.ipac',
            overwrite=True
        )
    else:
        tbl.write(
            '../output/tables/'+cat+'/'+str(n_regions)+'_distance_total.tbl',
            format='ascii.ipac',
            overwrite=True
        )
        


    skip = True
    if input('do you want to plot centroid distance vs YSO? (y/n) ').lower() == 'y':
        
        distribution_plots(n_regions,mode,cat)
    
if not skip:
    if input('do you want to plot centroid distance vs YSO? (y/n) ').lower() == 'y':
        
        print('(required file: (n_regions)_distance_(mode).tbl)' )
        n_regions = input('number of regions? (max 1017) ').lower()

        mode = '0'
        while (mode != 'match' and mode!= 'total'): 
            mode = input('choose sort (match/total): ').lower()

        distribution_plots(n_regions,mode,cat)


