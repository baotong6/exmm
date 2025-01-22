from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp  # 使用更快的插值方法

def crop_and_reproject_fits(input_fits, output_fits, row_start, row_end, col_start, col_end):
    """
    Crop a FITS file and reproject the cropped region with correct WCS coordinates.

    Parameters:
        input_fits (str): Path to the input FITS file.
        output_fits (str): Path to save the cropped FITS file.
        row_start (int): Starting row index (inclusive, 0-based).
        row_end (int): Ending row index (exclusive, 0-based).
        col_start (int): Starting column index (inclusive, 0-based).
        col_end (int): Ending column index (exclusive, 0-based).
    """
    # Load the input FITS file
    with fits.open(input_fits) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

        # Crop the data
        cropped_data = data[row_start:row_end, col_start:col_end]

        # Update the WCS for the cropped region
        new_wcs = WCS(naxis=2)
        new_wcs.wcs.crpix = [
            wcs.wcs.crpix[0] - col_start,
            wcs.wcs.crpix[1] - row_start,
        ]
        new_wcs.wcs.crval = wcs.wcs.crval
        new_wcs.wcs.cdelt = wcs.wcs.cdelt
        new_wcs.wcs.ctype = wcs.wcs.ctype
        new_wcs.wcs.cunit = wcs.wcs.cunit

        # Prepare target header
        target_header = new_wcs.to_header()
        shape_out = (row_end - row_start, col_end - col_start)
        print('123')
        # Reproject using interpolation for efficiency
        reprojected_data, footprint = reproject_interp((cropped_data, new_wcs), target_header, shape_out=shape_out)
        print('456')
        # Save to a new FITS file
        hdu = fits.PrimaryHDU(data=reprojected_data, header=target_header)
        hdu.writeto(output_fits, overwrite=True)

# Example usage
path = '/Users/baotong/data_GalDisc/data/'
input_fits = path + "GalDisc_ima_14_subdiv_smooth.fits.gz"
output_fits = path + "Galcen_14_subdiv_smooth_corrected.fits.gz"
row_start, row_end = 7000, 10000  # Row range
col_start, col_end = 7500, 10500  # Column range

crop_and_reproject_fits(input_fits, output_fits, row_start, row_end, col_start, col_end)
