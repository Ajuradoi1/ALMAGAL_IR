from astropy.table import Table
from YSO_statistics import YSO_plots_compare, region_plot
import numpy as np

#ALMA = input("Do you want to include ALMA data? (y/n): ").strip().lower() == "y" 

#cat = input("\ninput catalog name: ")

ALMA = True
cat = 'glimpse_s07'

if ALMA:
    spitzer_sources_alma = Table.read("../output/tables/"+cat+"/spitzer_sources_YSO_PLOTS_ALMA.tbl", format="ascii.ipac")
    alma_regions_alma = Table.read('../output/tables/'+cat+'/alma_regions_YSO_PLOTS_ALMA.tbl', format='ascii.ipac')
    alma_sources_alma = Table.read('../output/tables/'+cat+'/alma_sources_YSO_PLOTS_ALMA.tbl', format='ascii.ipac')


spitzer_sources = Table.read("../output/tables/"+cat+"/spitzer_sources_YSO_PLOTS.tbl", format="ascii.ipac")
alma_regions = Table.read('../output/tables/'+cat+'/alma_regions_YSO_PLOTS.tbl', format='ascii.ipac')
alma_sources = Table.read('../output/tables/'+cat+'/alma_sources_YSO_PLOTS.tbl', format='ascii.ipac')

#YSO_plots_compare(spitzer_sources,alma_regions,alma_sources,spitzer_sources_alma,alma_regions_alma,alma_sources_alma,ALMA,cat)

for i in range(5):
    region = np.random.choice(alma_regions['INTERNAL_ID'])

    region_plot(spitzer_sources,spitzer_sources_alma,cat,region)