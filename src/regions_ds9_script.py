# Loading libraries

import os
import pandas as pd
from astropy.table import Table

#catalog = 'glimpse2_v2cat'
#catalog = 'glimpse_s07'

def regions_ds9(alma_regions_file:str,alma_sources_file:str,catalog:str):

    alma_regions = Table.read("../alma_input/"+alma_regions_file, format="ascii.ipac")
    print(alma_regions[:2])
    print("\n\n\n")

    spitzer_sources = Table.read("../output/tables/spitzer_sources_"+catalog+".tbl", format="ascii.ipac")
    print(spitzer_sources[:2])
    print("\n\n\n")

    alma_sources = pd.read_csv('../alma_input/'+alma_sources_file, comment='\\', skiprows=[2], delimiter='\s+')
    print(alma_sources[6001:6004])
    print("\n\n\n")

    os.makedirs('../output/alma_sources', exist_ok=True)

    for region in alma_regions:

        alma_sources_cut = alma_sources[alma_sources['CLUMP'] == region['INTERNAL_ID']].reset_index(drop=True) # is a pandas object
        spitzer_sources_cut = spitzer_sources[spitzer_sources['internal_id_01'] == region['INTERNAL_ID']] # is a Table object

        os.makedirs(f'../output/alma_sources/{region["INTERNAL_ID"]}', exist_ok=True)
        
        fwhm_circ = 0.5

        if len(alma_sources_cut) > 0:    
            with open(f'../output/alma_sources/{region["INTERNAL_ID"]}/{region["INTERNAL_ID"]}_alma_sources.reg', 'w') as f:
                f.write('# Region file format: DS9 version 4.1\n')
                f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                f.write('icrs\n')
                
                for alma_index in range(len(alma_sources_cut)):
                    f.write('circle('+str(alma_sources_cut['RA'][alma_index])+','+str(alma_sources_cut['DEC'][alma_index])+','+str(fwhm_circ/3600)+') # text={'+str(alma_index+1)+'}\n')


        if len(spitzer_sources_cut) > 0:
            with open(f'../output/alma_sources/{region["INTERNAL_ID"]}/{region["INTERNAL_ID"]}_{catalog}_spitzer_sources.reg', 'w') as f:
                f.write('# Region file format: DS9 version 4.1\n')
                f.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                f.write('icrs\n')

                for spitzer_index in range(len(spitzer_sources_cut)):
                    f.write('circle('+str(spitzer_sources_cut['ra'][spitzer_index])+','+str(spitzer_sources_cut['dec'][spitzer_index])+','+str(fwhm_circ/3600)+')\n')
                    