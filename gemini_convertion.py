import numpy as np
from astropy.table import Table
from astropy.io import fits
import os

# --- Configuration ---
# Make sure these paths are correct for your environment
HEAD_FITS_PATH = "/global/cfs/cdirs/lsst/www/jolteon/data/FINAL2/JOLTEON_FINAL_0000_HEAD.FITS"
PHOT_FITS_PATH = "/global/cfs/cdirs/lsst/www/jolteon/data/FINAL2/JOLTEON_FINAL_0000_PHOT.FITS"

# Directory where individual SNANA FITS files will be saved
OUTPUT_SNANA_DIR = "snana_data"
# Name of the .LIST file to create
OUTPUT_LIST_FILE = "snana_data.LIST"

# --- Create output directory if it doesn't exist ---
os.makedirs(OUTPUT_SNANA_DIR, exist_ok=True)

# --- Load the full HEAD and PHOT tables ---
print(f"Loading HEAD data from: {HEAD_FITS_PATH}")
jolteon_head = Table(fits.open(HEAD_FITS_PATH)[1].data)
print(f"Loading PHOT data from: {PHOT_FITS_PATH}")
jolteon_phot = Table(fits.open(PHOT_FITS_PATH)[1].data)

print(f"Total supernovae in HEAD file: {len(jolteon_head)}")
print(f"Total photometric points in PHOT file: {len(jolteon_phot)}")

# --- Prepare a list to store paths of generated SNANA files ---
snana_file_paths = []

# --- Iterate through each supernova in the HEAD table ---
print("Processing individual supernovae...")
for idx, sn_head_row in enumerate(jolteon_head):
    snid = sn_head_row['SNID']
    # Extract light curve for the current SNID
    ptr_min = sn_head_row['PTROBS_MIN']
    ptr_max = sn_head_row['PTROBS_MAX']
    lc = jolteon_phot[ptr_min-1:ptr_max] # -1 because PTROBS_MIN/MAX are 1-indexed

    # --- Create the SNANA-style FITS file for this supernova ---
    # SNANA FITS files typically have two HDUs:
    # HDU1: Header info (like your HEAD data)
    # HDU2: Photometry data (like your PHOT data)

    # 1. Create a header table (HDU1)
    # Include PTROBS_MIN and PTROBS_MAX as BayeSN (via sncosmo) expects them.
    head_cols_to_keep = [
        'SNID', 'RA', 'DEC', 'MWEBV', 'MWEBV_ERR', 'SPECZ', 'SPECZ_ERR',
        'PHOTOZ', 'PHOTOZ_ERR', 'NOBS', 'PTROBS_MIN', 'PTROBS_MAX'
    ]
    # For a single row table, wrap each scalar value in a list
    snana_head_table = Table([[sn_head_row[col]] for col in head_cols_to_keep], names=head_cols_to_keep)

    # SNANA expects specific column names for photometry.
    # We will now create a 'BAND' column and ensure its contents are byte strings.
    band_mapping = {
        'g': 'LSST_g',
        'r': 'LSST_r',
        'i': 'LSST_i',
        'z': 'LSST_z',
        'y': 'LSST_y'
    }

    # Create the photometry table (HDU2)
    snana_phot_data = Table()
    snana_phot_data['MJD'] = lc['MJD']
    
    # --- IMPORTANT CHANGE HERE: Encoding band names to byte strings ---
    # Map original band names to LSST_g/r/i/z/y format, then encode to bytes
    encoded_bands = [band_mapping.get(b, b).encode('utf-8') for b in lc['BAND']]
    snana_phot_data['BAND'] = np.array(encoded_bands) # Ensure it's a NumPy array of byte strings
    
    snana_phot_data['FLUXCAL'] = lc['FLUXCAL']
    snana_phot_data['FLUXERR'] = lc['FLUXCALERR']
    snana_phot_data['KCOR'] = np.full(len(lc), np.nan) # Placeholder, fill with actual KCOR if you have it
    snana_phot_data['ZEROPT'] = lc['ZEROPT']

    # --- Construct the FITS HDUs ---
    primary_hdu = fits.PrimaryHDU()
    head_hdu = fits.BinTableHDU(snana_head_table, name='HEAD') # Name it 'HEAD'
    phot_hdu = fits.BinTableHDU(snana_phot_data, name='PHOT') # Name it 'PHOT'

    hdul = fits.HDUList([primary_hdu, head_hdu, phot_hdu])

    # Define output filename
    snana_filename = os.path.join(OUTPUT_SNANA_DIR, f"SN_{snid:06d}.FITS") # Use 6 digits for consistency
    hdul.writeto(snana_filename, overwrite=True) # overwrite=True for development

    snana_file_paths.append(os.path.abspath(snana_filename)) # Store absolute path

    if (idx + 1) % 100 == 0:
        print(f"  Processed {idx + 1} supernovae...")

print("\nAll supernovae processed.")

# --- Create the .LIST file ---
list_file_path = os.path.join(OUTPUT_SNANA_DIR, OUTPUT_LIST_FILE)
with open(list_file_path, 'w') as f:
    for path in snana_file_paths:
        f.write(f"{path}\n")

print(f"Created .LIST file at: {list_file_path}")
print(f"You can now point BayeSN's 'private_data_path' in your YAML to: {os.path.abspath(list_file_path)}")