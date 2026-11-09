"""Generate an annular ARF from a point-source ARF using CIAO CALDB REEF."""

import argparse
from pathlib import Path

from astropy.io import fits
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_arf", type=Path)
    parser.add_argument("output_arf", type=Path)
    parser.add_argument("--theta-arcmin", type=float, required=True)
    parser.add_argument("--phi-deg", type=float, required=True)
    parser.add_argument("--reef", type=Path, help="Explicit CALDB REEF file")
    args = parser.parse_args()
    if not np.isfinite(args.theta_arcmin) or args.theta_arcmin < 0 or not np.isfinite(args.phi_deg):
        parser.error("Supply finite nominal off-axis angles from this observation")
    csv_path = args.output_arf.with_name("reef_fraction.csv")
    if args.output_arf.exists() or csv_path.exists():
        parser.error("Use a new output directory to preserve existing ARF and REEF products")

    import psf
    if args.reef is None:
        import caldb4
        candidates = {Path(entry.split("[", 1)[0]) for entry in
                      caldb4.Caldb(telescope="CHANDRA", product="REEF").search}
        if len(candidates) != 1:
            parser.error("CALDB REEF query was ambiguous. Supply the intended file with --reef")
        reef, = candidates
    else:
        reef = args.reef

    with fits.open(args.input_arf, memmap=False) as hdus:
        original = hdus["SPECRESP"]
        if "PSF_FRAC" in original.columns.names:
            raise ValueError("Input ARF already contains an aperture fraction")
        energy = (original.data["ENERG_LO"] + original.data["ENERG_HI"]) / 2
        calibration = psf.psfInit(str(reef))
        try:
            fraction = np.array([
                psf.psfFrac(calibration, float(e), args.theta_arcmin, args.phi_deg, 2.0)
                - psf.psfFrac(calibration, float(e), args.theta_arcmin, args.phi_deg, 0.75)
                for e in energy
            ])
        finally:
            psf.psfClose(calibration)
        if not np.all(np.isfinite(fraction) & (fraction > 0) & (fraction <= 1)):
            raise ValueError("Nonphysical REEF annular fractions")
        column = fits.Column(name="PSF_FRAC", format="E", array=fraction)
        corrected = fits.BinTableHDU.from_columns(original.columns + fits.ColDefs([column]),
                                                 header=original.header, name=original.name)
        corrected.data["SPECRESP"] *= fraction
        corrected.header.add_history("CALDB REEF ECF(2.0 arcsec) minus ECF(0.75 arcsec)")
        corrected.header.add_history(f"REEF file: {reef.name}")
        corrected.header.add_history(f"theta={args.theta_arcmin} arcmin, phi={args.phi_deg} deg")
        hdus[hdus.index_of("SPECRESP")] = corrected
        hdus.writeto(args.output_arf, overwrite=False, checksum=True)
    np.savetxt(csv_path, np.column_stack([energy, fraction]), delimiter=",", comments="",
               header="energy_keV,annular_fraction")
    print(f"f_PSF(4 keV) = {np.interp(4.0, energy, fraction):.8f}")


if __name__ == "__main__":
    main()
