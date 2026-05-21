from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np


def twomass_to_spitzer(spitzer_sources, cat):

    twomass = Table.read("../output/tables/spitzer_2MASS_fp_psc.tbl", format="ipac")

    twomass_coords = SkyCoord(twomass["ra"], twomass["dec"], unit="deg", frame="icrs")

    spitzer_coords = SkyCoord(
        spitzer_sources["ra"], spitzer_sources["dec"], unit="deg", frame="icrs"
    )

    idx, sep, _ = spitzer_coords.match_to_catalog_sky(twomass_coords)

    radius = 1.0 * u.arcsec
    mask = sep < radius

    spitzer_sources["has_twomass"] = mask

    for col in ["j_m", "j_msigcom", "h_m", "h_msigcom", "k_m", "k_msigcom"]:
        spitzer_sources[col] = np.nan
        spitzer_sources[col][mask] = np.array(twomass[col][idx[mask]], dtype=float)
        # spitzer_sources[col].dtype

    spitzer_sources.write(
        "../output/tables/spitzer_sources_" + cat + "_2mass.tbl",
        format="ascii.ipac",
        overwrite=True,
    )

    return spitzer_sources
