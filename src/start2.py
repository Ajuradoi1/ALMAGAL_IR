from astropy.table import Table
from YSO_statistics import YSO_plots


ALMA = input("Do you want to include ALMA data? (y/n): ").strip().lower() == "y" 

if ALMA:
    spitzer_sources = Table.read("../output/tables/spitzer_sources_YSO_PLOTS_ALMA.tbl", format="ascii.ipac")
    alma_regions = Table.read('../output/tables/alma_regions_YSO_PLOTS_ALMA.tbl', format='ascii.ipac')
    alma_sources = Table.read('../output/tables/alma_sources_YSO_PLOTS_ALMA.tbl', format='ascii.ipac')
else:
    spitzer_sources = Table.read("../output/tables/spitzer_sources_YSO_PLOTS.tbl", format="ascii.ipac")
    alma_regions = Table.read('../output/tables/alma_regions_YSO_PLOTS.tbl', format='ascii.ipac')
    alma_sources = Table.read('../output/tables/alma_sources_YSO_PLOTS.tbl', format='ascii.ipac')

cat = 'glimpse_s07_ALMA' if ALMA else 'glimpse_s07'

YSO_plots(spitzer_sources=spitzer_sources, alma_regions=alma_regions, alma_sources=alma_sources, alma_proceeding=ALMA, cat=cat)