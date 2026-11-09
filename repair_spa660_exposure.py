"""Legacy exposure-only check for products_fast; use reextract_spa660.py for channel validation."""

import json
from pathlib import Path

import numpy as np
from astropy.io import fits


def matching_exposure(value, exposure):
    return np.isclose(value, exposure, rtol=0.01, atol=2)


def main():
    folder = Path("products_fast")
    event = Path(Path("products/event_path.txt").read_text().strip())
    with fits.open(event) as hdus:
        events = hdus["EVENTS"]
        if int(events.header["OBS_ID"]) != 14946 or set(events.data["CCD_ID"]) != {7}:
            raise ValueError("Expected single-chip CCD 7 events from ObsID 14946")
        gti, = (
            hdu.data for hdu in hdus[2:]
            if hdu.header.get("HDUNAME") == "GTI7"
            or (hdu.name == "GTI" and hdu.header.get("CCD_ID", hdu.ver) == 7)
        )
        start, stop = gti["START"], gti["STOP"]
        dtcor = float(events.header["DTCOR"])
        if not len(start) or not 0 < dtcor <= 1 or np.any(stop <= start) or np.any(start[1:] < stop[:-1]):
            raise ValueError("Invalid CCD 7 good-time intervals or live-time factor")
        exposure = float(np.sum(stop - start) * dtcor)
        header_live = float(events.header.get("LIVETIME7", events.header.get("LIVETIME", 0)))

    aspect = fits.getdata("products/asphist_ccd7.fits", "ASPHIST")
    aspect_live = float(np.sum(aspect["DURATION"]))
    print(f"CCD 7 live time: GTI × DTCOR = {exposure:.3f} s, "
          f"event header = {header_live:.3f} s, aspect histogram = {aspect_live:.3f} s")
    if not np.isfinite(exposure) or exposure <= 0 or not matching_exposure(aspect_live, exposure):
        raise ValueError("The aspect exposure does not match the CCD 7 GTI-derived exposure")
    if header_live > 0 and not matching_exposure(header_live, exposure):
        raise ValueError("The event header live time does not match the CCD 7 GTIs")

    bounds = fits.getdata(folder / "magnetar.rmf", "EBOUNDS")
    energies = (bounds["E_MIN"] + bounds["E_MAX"]) / 2
    channels = bounds["CHANNEL"][(energies >= 1) & (energies < 10)]
    measured = json.loads((folder / "event_check/event_summary.json").read_text())
    records = {}
    for name, expected in (
        ("magnetar.pi", measured["source_counts_1_10"]),
        ("magnetar_grp.pi", measured["source_counts_1_10"]),
        ("magnetar_bkg.pi", measured["background_counts_1_10"]),
    ):
        path = folder / name
        with fits.open(path) as hdus:
            pha = hdus["SPECTRUM"]
            counts = pha.data["COUNTS"]
            in_band = np.isin(pha.data["CHANNEL"], channels)
            total = int(counts[in_band].sum())
            old = float(pha.header.get("EXPOSURE", 0))
            scale = float(pha.header["BACKSCAL"])
            print(f"{name}: EXPOSURE={old:g} s, BACKSCAL={scale:.6g}, "
                  f"1–10 keV PI counts={total}, event-energy counts={expected}")
            if total <= 0 or abs(total - expected) > max(25, 0.10 * expected):
                raise ValueError(f"{name} counts disagree with the independent event extraction")
            if not np.isfinite(scale) or scale <= 0:
                raise ValueError(f"{name} has invalid BACKSCAL")
            if not np.isfinite(old) or old < 0 or (old > 0 and not matching_exposure(old, exposure)):
                raise ValueError(f"{name} exposure conflicts with the verified live time")
            records[name] = (old, scale)

    ratio = records["magnetar_grp.pi"][1] / records["magnetar_bkg.pi"][1]
    if not np.isclose(ratio, measured["area_ratio"], rtol=0.15):
        raise ValueError(f"BACKSCAL ratio {ratio:g} disagrees with the masked aperture area")
    print(f"Verified source/background BACKSCAL ratio: {ratio:.5f}")
    for name, (old, _) in records.items():
        if old > 0:
            continue
        with fits.open(folder / name, mode="update") as hdus:
            pha = hdus["SPECTRUM"]
            pha.header["EXPOSURE"] = (exposure, "CCD 7 exposure in seconds")
            pha.header["LIVETIME"] = exposure
            if "COUNT_RATE" in pha.columns.names:
                pha.data["COUNT_RATE"] = pha.data["COUNTS"] / exposure
            pha.header.add_history(f"Original EXPOSURE was {old:g} s.")
            pha.header.add_history("Exposure restored from CCD 7 GTIs times DTCOR.")
            pha.header.add_history("Checked against the CCD 7 aspect histogram duration.")
            pha.add_checksum()
    print(f"All three spectra have positive exposure consistent with {exposure:.3f} s.")
    print("This legacy repair does not verify individual PI channels.")


if __name__ == "__main__":
    main()
