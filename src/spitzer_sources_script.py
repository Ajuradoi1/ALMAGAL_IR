# Loading libraries

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import requests
import os


def spitzer_sources(alma_regions_file:str,catalog:str):

    # filename = "almaprops_2024_all_20250219.tab"
    filename = alma_regions_file
    alma_regions = Table.read(f"../alma_input/{filename}", format="ascii.ipac")
    print(alma_regions)
    print("\ngenerating coordinates file...\n")

    os.makedirs("../output", exist_ok=True)
    os.makedirs("../output/tables",exist_ok=True)

    with open("../output/tables/alma_regions_coordinates.tbl", "w", newline="\n") as f:
        f.write("|INTERNAL_ID         |ra         |dec       |\n")
        f.write("|char                |float      |float     |\n")
        for row in alma_regions:
            coords = SkyCoord(
                row["GLON"], row["GLAT"], unit=(u.deg, u.deg), frame="galactic"
            )
            f.write(
                f"{row['INTERNAL_ID']:>20}  {coords.icrs.ra.deg:>3.4f}    {coords.icrs.dec.deg:>3.4f}\n"
            )

    # "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=glimpse_s07&spatial=upload&radius=60"

    print("querying Spitzer data...\n")

    url = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
    file = {"filename": open("../output/tables/alma_regions_coordinates.tbl", "rb")}

    #catalog = 'glimpse2_v2cat'
    #catalog = 'glimpse_s07'

    data = {
        "catalog": catalog,
        "spatial": "Upload",
        "uradius": 60,
        "outfmt": 1
        }

    response = requests.post(url, files=file, data=data)
    with open(f"../output/tables/spitzer_sources_{catalog}.tbl", "w") as f:
        f.write(response.text)

    print("done!\n")

