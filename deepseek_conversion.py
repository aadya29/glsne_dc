import astropy.table as at
from astropy.io import fits
import os
import numpy as np
import sys

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Configuration - use paths relative to script location
output_dir = os.path.join(script_dir, "JOLTEON_SNANA_ASCII")
list_filename = os.path.join(script_dir, "JOLTEON_SNANA_ASCII.LIST")

print(f"Script directory: {script_dir}")
print(f"Output directory: {output_dir}")
print(f"List file: {list_filename}")

# Create output directory with error handling
try:
    os.makedirs(output_dir, exist_ok=True)
    print(f"Created output directory: {output_dir}")
except Exception as e:
    print(f"Error creating directory {output_dir}: {e}")
    sys.exit(1)

# Verify write permissions
if not os.access(output_dir, os.W_OK):
    print(f"Error: No write permission for directory {output_dir}")
    sys.exit(1)

# Load data
head_path = "/global/cfs/cdirs/lsst/www/jolteon/data/FINAL2/JOLTEON_FINAL_0000_HEAD.FITS"
phot_path = "/global/cfs/cdirs/lsst/www/jolteon/data/FINAL2/JOLTEON_FINAL_0000_PHOT.FITS"

try:
    head = at.Table(fits.open(head_path)[1].data)
    phot = at.Table(fits.open(phot_path)[1].data)
    print(f"Loaded {len(head)} supernovae from header table")
    print(f"Loaded {len(phot)} photometry points")
except Exception as e:
    print(f"Error loading FITS files: {e}")
    sys.exit(1)

# Create list file
try:
    with open(list_filename, 'w') as list_f:
        for i, sn in enumerate(head):
            snid = sn['SNID']
            ptr_min = sn['PTROBS_MIN']
            ptr_max = sn['PTROBS_MAX']
            
            # Handle potential index issues
            if ptr_min < 1 or ptr_max > len(phot):
                print(f"Warning: Invalid PTROBS range for SNID {snid} ({ptr_min}-{ptr_max})")
                continue
                
            lc = phot[ptr_min-1:ptr_max]
            
            lc_filename = f"{snid}.DAT"
            lc_filepath = os.path.join(output_dir, lc_filename)
            list_f.write(f"{lc_filename}\n")
            
            with open(lc_filepath, 'w') as f:
                # Write SNANA-style header
                f.write(f"# ================================\n")
                f.write(f"# SNID: {snid}\n")
                f.write(f"# RA: {sn['RA']}\n")
                f.write(f"# DEC: {sn['DEC']}\n")
                f.write(f"# ================================\n")
                f.write(f"NVAR: 4\n")
                f.write(f"VARNAMES: MJD FLT MAG MAGERR\n")
                
                for obs in lc:
                    # Convert filter name
                    flt = obs['BAND'] + '_LSST'
                    
                    # Convert flux to magnitude
                    flux = obs['FLUXCAL']
                    fluxerr = obs['FLUXCALERR']
                    
                    # Skip invalid flux values
                    if flux <= 0:
                        continue
                    
                    # Calculate magnitude and error
                    mag = obs['ZEROPT'] - 2.5 * np.log10(flux)
                    magerr = (2.5 / np.log(10)) * (fluxerr / flux)
                    
                    f.write(f"{obs['MJD']:.5f}  {flt}  {mag:.5f}  {magerr:.5f}\n")
            
            # Print progress
            if (i + 1) % 100 == 0:
                print(f"Processed {i+1}/{len(head)} supernovae")

    print(f"\nSuccessfully created {len(head)} light curves")
    print(f"List file: {list_filename}")
    print(f"Light curve directory: {output_dir}")
    print("\nImportant: Make sure your BayeSN configuration includes:")
    print("version_photometry: JOLTEON_SNANA_ASCII")

except Exception as e:
    print(f"Error during processing: {e}")
    sys.exit(1)