"""Select ObsID 14946 events and verify CIAO PI spectra before fitting."""

import argparse
import json
from pathlib import Path
import shutil
import warnings

from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
import numpy as np

from inspect_events import MAGNETAR, MASKS, aperture_area


OUTPUT = Path("products_verified")
FAST = Path("products_fast")
APERTURES = {"source": (0.75, 2.0), "background": (2.5, 8.0)}


def select_events(output, select_only):
    output.mkdir(exist_ok=True)
    evt = Path(Path("products/event_path.txt").read_text().strip())
    with fits.open(evt) as hdus:
        events = hdus["EVENTS"]
        header = events.header
        data = events.data
        if int(header["OBS_ID"]) != 14946 or set(data["ccd_id"]) != {7}:
            raise ValueError("Expected CCD 7 events from ObsID 14946")
        gti, = (
            hdu.data for hdu in hdus[2:]
            if hdu.header.get("HDUNAME") == "GTI7"
            or (hdu.name == "GTI" and hdu.header.get("CCD_ID", hdu.ver) == 7)
        )
        start, stop = gti["START"], gti["STOP"]
        dtcor = float(header["DTCOR"])
        if not len(start) or not 0 < dtcor <= 1 or np.any(stop <= start) or np.any(start[1:] < stop[:-1]):
            raise ValueError("Invalid CCD 7 good-time intervals or live-time factor")
        ontime = float(np.sum(stop - start))
        exposure = ontime * dtcor
        hist = fits.getdata("products/asphist_ccd7.fits", "ASPHIST")
        aspect_exposure = float(np.sum(hist["DURATION"]))
        if not np.isfinite(exposure) or exposure <= 0 or not np.isclose(aspect_exposure, exposure, rtol=0.01):
            raise ValueError("CCD 7 good-time exposure disagrees with the aspect histogram")

        column_names = [name.lower() for name in events.columns.names]
        xy_cols = [column_names.index(axis) + 1 for axis in ("x", "y")]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(header, keysel=["pixel"], colsel=xy_cols, naxis=2)
        scale = np.abs(wcs.wcs.cdelt) * 3600
        if not np.isclose(scale[0], scale[1], rtol=1e-8):
            raise ValueError("Circular SKY regions require square pixels")
        origin = wcs.all_world2pix([MAGNETAR], 1)[0]
        x = (data["x"] - origin[0]) * scale[0]
        y = (data["y"] - origin[1]) * scale[1]
        radius = np.hypot(x, y)
        offset_masks = [
            (*((wcs.all_world2pix([[ra, dec]], 1)[0] - origin) * scale), size)
            for ra, dec, size in MASKS
        ]
        areas = {
            "source": np.pi * (2.0**2 - 0.75**2),
            "background": aperture_area(offset_masks),
        }
        meta = {
            "obsid": 14946, "event_file": str(evt), "exposure_s": exposure,
            "aspect_exposure_s": aspect_exposure, "ontime_s": ontime, "dtcor": dtcor,
            "areas_arcsec2": areas, "backscale": {},
        }
        for label, (inner, outer) in APERTURES.items():
            selection = (radius >= inner) & (radius < outer)
            region = [f"annulus({origin[0]:.12f},{origin[1]:.12f},"
                      f"{inner/scale[0]:.12f},{outer/scale[0]:.12f})"]
            if label == "background":
                for mx, my, size in offset_masks:
                    selection &= np.hypot(x - mx, y - my) >= size
                    region.append(f"-circle({origin[0]+mx/scale[0]:.12f},"
                                  f"{origin[1]+my/scale[1]:.12f},{size/scale[0]:.12f})")
            selected = data[selection]
            energy = selected["energy"]
            band_count = int(np.count_nonzero((energy >= 1000) & (energy < 10000)))
            if not band_count:
                raise ValueError(f"The {label} event selection is empty")
            meta[f"{label}_events_1_10"] = band_count
            meta["backscale"][label] = areas[label] / (scale.prod() * 8192**2)
            events.data = selected
            hdus.writeto(output / f"{label}_events.fits", overwrite=False, checksum=True)
            (output / f"{label}_physical.reg").write_text(
                "physical\n" + "\n".join(region) + "\n", encoding="ascii", newline="\n"
            )
            print(f"{label}: {len(selected)} total events, {band_count} in 1–10 keV")
    if not select_only:
        independent = json.loads((FAST / "event_check/event_summary.json").read_text())
        for label in APERTURES:
            if meta[f"{label}_events_1_10"] != independent[f"{label}_counts_1_10"]:
                raise ValueError(f"{label} selection differs from the independent event check")
        for filename in ("magnetar.arf", "magnetar.corr.arf", "magnetar.rmf", "reef_fraction.csv"):
            shutil.copyfile(FAST / filename, output / filename)
        shutil.copytree(FAST / "event_check", output / "event_check", dirs_exist_ok=True)
    (output / "extraction_summary.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"Verified CCD 7 exposure: {exposure:.3f} s")


def verify_pi(output):
    meta = json.loads((output / "extraction_summary.json").read_text())
    if meta["obsid"] != 14946 or not 0 < meta["exposure_s"] < np.inf:
        raise ValueError("Invalid ObsID or exposure in extraction summary")
    for label, filename in (("source", "magnetar.pi"), ("background", "magnetar_bkg.pi")):
        pi = fits.getdata(output / f"{label}_events.fits", "EVENTS", memmap=False)["PI"]
        expected = np.bincount(pi.astype(int))
        with fits.open(output / filename, mode="update") as hdus:
            spectrum = hdus["SPECTRUM"]
            channel = spectrum.data["CHANNEL"].astype(int)
            observed = spectrum.data["COUNTS"]
            in_range = (channel >= 0) & (channel < len(expected))
            selected_counts = np.zeros_like(channel)
            selected_counts[in_range] = expected[channel[in_range]]
            if not np.array_equal(selected_counts, observed) or observed.sum() <= 0:
                raise ValueError(f"{filename} PI counts disagree with the selected events")
            if observed.sum() != len(pi):
                raise ValueError(f"{filename} does not cover every selected PI channel")
            header = spectrum.header
            old_exposure = float(header.get("EXPOSURE", 0))
            if not np.isfinite(old_exposure) or old_exposure < 0:
                raise ValueError(f"{filename} has invalid exposure")
            if old_exposure and not np.isclose(old_exposure, meta["exposure_s"], rtol=0.01):
                raise ValueError(f"{filename} has a conflicting positive exposure")
            header["EXPOSURE"] = (meta["exposure_s"], "CCD 7 GTI live time in seconds")
            header["LIVETIME"] = meta["exposure_s"]
            header["ONTIME"] = meta["ontime_s"]
            header["BACKSCAL"] = (meta["backscale"][label], "Aperture area / ACIS sky area")
            header["AREASCAL"] = 1.0
            header["POISSERR"] = True
            if "COUNT_RATE" in spectrum.columns.names:
                spectrum.data["COUNT_RATE"] = observed / meta["exposure_s"]
            if label == "source":
                header["BACKFILE"] = "magnetar_bkg.pi"
                header["RESPFILE"] = "magnetar.rmf"
                header["ANCRFILE"] = "magnetar.corr.arf"
            header.add_history("CIAO dmextract run on a validated physical event selection.")
            header.add_history("Exposure from CCD 7 GTIs and DTCOR, checked with ASPHIST.")
            header.add_history("BACKSCAL from geometric aperture area including source masks.")
            spectrum.add_checksum()
            print(f"Verified {filename}: {int(observed.sum())} PI counts, "
                  f"EXPOSURE={meta['exposure_s']:.3f} s, BACKSCAL={meta['backscale'][label]:.6g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("select", "verify"))
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--select-only", action="store_true",
                        help="Write selected events and metadata without the legacy products_fast inputs")
    args = parser.parse_args()
    if args.stage == "select":
        select_events(args.output, args.select_only)
    else:
        if args.select_only:
            parser.error("--select-only is only supported by the select stage")
        verify_pi(args.output)
