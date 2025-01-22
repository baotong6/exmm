import csv
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from tqdm import tqdm

def estimate_false_match_rate(ra, dec, radius_arcsec, large_area_radius_deg=0.05):
    """
    Estimate the false match rate for matching X-ray sources with Gaia sources.

    Parameters:
        ra (float): Right Ascension (RA) of the center point (degrees).
        dec (float): Declination (Dec) of the center point (degrees).
        radius_arcsec (float): Matching radius in arcseconds.
        large_area_radius_deg (float): Radius of the larger area to estimate Gaia density (degrees).

    Returns:
        tuple: (source_density_per_sq_sec, false_match_rate)
    """
    radius_deg = radius_arcsec / 3600.0
    coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

    try:
        result = Vizier(row_limit=-1, columns=["*"]).query_region(
            coord, radius=large_area_radius_deg * u.degree, catalog="I/355"
        )
    except Exception as e:
        print(f"Error during Vizier query: {e}")
        return None, None

    if len(result) == 0:
        print("No Gaia sources found in the large search region.")
        return None, None

    gaia_table = result[0]
    number_of_sources = len(gaia_table)
    area_of_large_search = np.pi * (large_area_radius_deg ** 2)

    source_density_per_sq_deg = number_of_sources / area_of_large_search
    source_density_per_sq_sec = source_density_per_sq_deg / 3600 ** 2
    area_of_matching_radius = np.pi * (radius_deg ** 2)
    expected_sources_in_radius = source_density_per_sq_deg * area_of_matching_radius
    false_match_rate = 1 - np.exp(-expected_sources_in_radius)

    return source_density_per_sq_sec, false_match_rate

def matches_togaia_with_selected_columns(srclist, output_file):
    """
    Match ESASS coordinates with Gaia catalog and append match results to selected columns from original data.
    Each source is queried separately to limit the region queried.

    Parameters:
        srclist (pd.DataFrame): DataFrame containing the source list with selected columns.
        output_file (str): Path to save the output CSV file.
    """
    selected_columns = [
        'e_id', 'e_ra', 'e_dec', 'xmm_index', 'xmm_ra', 'xmm_dec',
        'separation_arcsec', 'hamstarindex', 'Fx','CTP_ID', 'CTP_SEP'
    ]
    srclist = srclist[selected_columns]
    srclist['CTP_ID'] = srclist['CTP_ID'].astype(str)  # Ensure CTP_ID is a string

    ra = srclist.apply(
        lambda row: row["xmm_ra"] if row["separation_arcsec"] <= 17 else row["e_ra"], axis=1
    ).values
    dec = srclist.apply(
        lambda row: row["xmm_dec"] if row["separation_arcsec"] <= 17 else row["e_dec"], axis=1
    ).values

    coord_e = SkyCoord(ra * u.degree, dec * u.degree)
    distance_limit = 20.0 * u.arcsec

    all_rows = []
    with open(output_file, mode='w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        header = selected_columns + [
            "gaia_index", "gaia_ra", "gaia_dec", "gaia_separation_arcsec",
            "gaia_parallax", "gaia_parallax_error", "gaia_pm", "gaia_pmra_error",
            "gaia_pmdec_error", "gaia_gmag", "gaia_bpmag", "gaia_rpmag",
            "source_density_per_sq_sec", "false_match_rate"
        ]
        csvwriter.writerow(header)

        print("Matching ESASS sources with Gaia catalog...")
        for i, source in tqdm(enumerate(coord_e), total=len(coord_e)):
            coord_e_sky = SkyCoord(ra[i] * u.degree, dec[i] * u.degree)
            try:
                result = Vizier(
                    row_limit=-1,
                    columns=["*"]
                ).query_region(coord_e_sky, radius=0.01 * u.degree, catalog="I/355")
            except Exception as e:
                print(f"Error during Vizier query for source {i}: {e}")
                continue

            density, false_match_rate = estimate_false_match_rate(ra[i], dec[i], 20.0)
            if len(result) == 0:
                row = srclist.iloc[i].tolist() + [-1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, density, false_match_rate]
                csvwriter.writerow(row)
                continue

            gaia_table = result[0]
            gaia_table['Source'] = gaia_table['Source'].astype(str)  # Ensure Gaia Source is a string
            coord_gaia = SkyCoord(gaia_table["RAJ2000"], gaia_table["DEJ2000"], unit="deg")
            separations = coord_e_sky.separation(coord_gaia)
            idxc = np.where(separations < distance_limit)[0]

            if len(idxc) > 0:
                for j in idxc:
                    row = srclist.iloc[i].tolist() + [
                        gaia_table["Source"][j],
                        coord_gaia[j].ra.deg,
                        coord_gaia[j].dec.deg,
                        separations[j].arcsecond,
                        gaia_table["Plx"][j],
                        gaia_table["e_Plx"][j],
                        gaia_table["PM"][j],
                        gaia_table["e_pmRA"][j],
                        gaia_table["e_pmDE"][j],
                        gaia_table["Gmag"][j],
                        gaia_table["BPmag"][j],
                        gaia_table["RPmag"][j],
                        density,
                        false_match_rate
                    ]
                    csvwriter.writerow(row)
                    all_rows.append(row)
            else:
                row = srclist.iloc[i].tolist() + [-1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, density, false_match_rate]
                csvwriter.writerow(row)
                all_rows.append(row)

    print(f"Results saved to {output_file}")
    excel_file = output_file[:-3] + 'xlsx'
    df = pd.DataFrame(all_rows, columns=header)
    df.to_excel(excel_file, index=False)  # Save as Excel file without row index
    print(f"Results also saved to {excel_file}")

if __name__ == '__main__':
    path = '/Users/baotong/data_GalDisc/data/match_e_xmm/'
    srclist = pd.read_excel(path + 'e_xmmdr14s_match_all_starinfo.xlsx')
    output_file = path + "matched_gaia_results_optimized_20sec.csv"
    matches_togaia_with_selected_columns(srclist, output_file)
