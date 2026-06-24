# Loading libraries

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import requests
import os


def spitzer_data(alma_regions_file: str, catalog: str, type: str):

    filename = alma_regions_file
    alma_regions = Table.read(f"../alma_input/{filename}", format="ascii.ipac")


    os.makedirs("../output", exist_ok=True)
    os.makedirs("../output/tables", exist_ok=True)
    

    if not os.path.isfile("../output/tables/alma_regions_coordinates.tbl"):
        print("\ngenerating coordinates file...\n")
        with open(
            "../output/tables/alma_regions_coordinates.tbl", "w", newline="\n"
        ) as f:
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

    print("querying spitzer data...\n")

    url = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
    file = {"filename": open("../output/tables/alma_regions_coordinates.tbl", "rb")}

    data = {"catalog": catalog, "spatial": "Upload", "uradius": 60, "outfmt": 1}

    response = requests.post(url, files=file, data=data)

    if (
        response.text.split(" list")[0]
        == '[struct stat="ERROR", msg="The catalog is not on the'
    ):
        print("invalid catalog name, please use a spitzer one.\n")
        exit()
    else:
        
        if catalog != 'fp_psc':
            
            os.makedirs("../output/tables/"+catalog, exist_ok=True)
            os.makedirs("../output/figures", exist_ok=True)
            os.makedirs("../output/figures/"+catalog, exist_ok=True)

            with open("../output/tables/"+catalog+f"/spitzer_{type}_{catalog}_download.tbl", "w") as f:
                f.write(response.text)

                print("done!\n")
        else:
            
            with open(f"../output/tables/spitzer_{type}_{catalog}_download.tbl", "w") as f:
                f.write(response.text)

                print("done!\n")
