"""Read-only audit of the saved ObsID 14946 products. Does not run CIAO or XSPEC."""

import argparse
import json
import math
from pathlib import Path

from astropy.io import fits
import numpy as np


def require(condition, description):
    if not condition:
        raise ValueError(description)


def read_table(path, extension):
    with fits.open(path, memmap=False) as hdus:
        hdus.verify("exception")
        for hdu in hdus:
            if hdu.header.get("CHECKSUM"):
                require(hdu.verify_checksum() == 1, f"Invalid checksum: {path.name}/{hdu.name}")
            datasum = str(hdu.header.get("DATASUM", "")).strip()
            if datasum:
                require(hdu.verify_datasum() == 1, f"Invalid data checksum: {path.name}/{hdu.name}")
        return hdus[extension].header.copy(), hdus[extension].data.copy()


def verify_pha(directory, name, events):
    header, spectrum = read_table(directory / name, "SPECTRUM")
    require(int(header["OBS_ID"]) == 14946, f"Incorrect ObsID: {name}")
    require(header["CHANTYPE"] == "PI" and header["POISSERR"], f"Invalid PI metadata: {name}")
    require(all(0 < header[key] < math.inf for key in ("EXPOSURE", "BACKSCAL", "AREASCAL")),
            f"Invalid exposure or scale: {name}")
    channel = spectrum["CHANNEL"].astype(int)
    require(np.all(channel >= 0) and np.all(np.diff(channel) == 1),
            f"Invalid channel grid: {name}")
    pi = events["PI"].astype(int)
    require(len(pi) and np.all((pi >= channel[0]) & (pi <= channel[-1])),
            f"Empty event list or PI outside PHA grid: {name}")
    counts = np.bincount(pi - channel[0], minlength=len(channel))
    require(np.array_equal(spectrum["COUNTS"], counts), f"Event/PHA mismatch: {name}")
    return header, spectrum


def verify_responses(directory):
    raw_header, unc = read_table(directory / "magnetar.arf", "SPECRESP")
    cor_header, corr = read_table(directory / "magnetar.corr.arf", "SPECRESP")
    require(int(raw_header["OBS_ID"]) == int(cor_header["OBS_ID"]) == 14946,
            "ARF observation does not match the event lists")
    require("PSF_FRAC" not in unc.names, "Uncorrected ARF contains a PSF_FRAC column")
    _, matrix = read_table(directory / "magnetar.rmf", "MATRIX")
    for key in ("ENERG_LO", "ENERG_HI"):
        require(np.array_equal(unc[key], corr[key]) and np.array_equal(unc[key], matrix[key]),
                "ARF/RMF grid mismatch")
    energy = (corr["ENERG_LO"] + corr["ENERG_HI"]) / 2
    fraction = corr["PSF_FRAC"]
    require(np.all(np.isfinite(fraction) & (fraction > 0) & (fraction <= 1)),
            "Nonphysical aperture fraction")
    require(np.allclose(unc["SPECRESP"] * fraction, corr["SPECRESP"], rtol=2e-6, atol=0),
            "Corrected ARF differs from the uncorrected ARF times PSF_FRAC")
    band = (energy >= 1) & (energy <= 10)
    require(np.any(band) and corr["ENERG_LO"][0] <= 1 and corr["ENERG_HI"][-1] >= 10,
            "Response does not cover 1–10 keV")
    reef = np.loadtxt(directory / "reef_fraction.csv", delimiter=",", skiprows=1)
    require(reef.shape == (len(energy), 2) and np.allclose(reef[:, 0], energy) and
            np.allclose(reef[:, 1], fraction), "REEF table differs from the corrected ARF")
    f4 = float(np.interp(4, energy, fraction))
    print(f"ARF: f(4 keV)={f4:.9f}, reciprocal={1/f4:.6f}, "
          f"1–10 keV range={fraction[band].min():.8f}–{fraction[band].max():.8f}")
    print("ARF multiplication verified numerically. Calibration generation NOT VERIFIED.")


def verify_fit_artifacts(directory):
    results = json.loads((directory / "fit_results.json").read_text(encoding="utf-8"))
    distance = float(results["distance_kpc"])
    factor = 4 * math.pi * (distance * 3.08567758149e21)**2
    expected = {f"tbabs_{model}_{column}" for model in ("bbody", "powerlaw")
                for column in ("free", "fixed")}
    expected.update({"phabs_bbody_free", "phabs_powerlaw_free"})
    require(expected == set(results["fits"]), "Missing model fit records")
    for name, fit in results["fits"].items():
        flux = fit["unabsorbed_flux_1_10"]
        lum = fit["luminosity_1_10"]
        require(math.isclose(lum["value"], factor * flux["value"], rel_tol=1e-10),
                f"Flux/luminosity mismatch: {name}")
        bounds = flux["interval_68"]
        require((bounds is None and lum["interval_68"] is None) or
                (bounds is not None and lum["interval_68"] is not None and
                 np.allclose(np.asarray(bounds) * factor, lum["interval_68"], rtol=1e-10, atol=0)),
                f"Luminosity interval mismatch: {name}")
        session = (directory / f"{name}.xcm").read_text(encoding="utf-8").splitlines()
        tokens = {tuple(line.split()) for line in session}
        for command, filename in (("data", "magnetar_grp.pi"), ("response", "magnetar.rmf"),
                                  ("arf", "magnetar.corr.arf"), ("backgrnd", "magnetar_bkg.pi")):
            require(any((command, index, filename) in tokens for index in ("1", "1:1")),
                    f"Incorrect saved {command} association: {name}")
        require(("statistic", "chi") in tokens and ("abund", "angr") in tokens and
                ("xsect", "bcmc") in tokens and ("xset", "TBABSVERSION", "2") in tokens,
                f"Unexpected XSPEC settings: {name}")
        for key in ("nh_1e22_cm2", "shape", "norm", "absorbed_flux_1_10", "unabsorbed_flux_1_10"):
            entry = fit[key]
            if entry["interval_68"] is not None:
                low, high = entry["interval_68"]
                require(entry.get("xspec_error_status") == "FFFFFFFFF" and
                        low < entry["value"] < high, f"Unusable saved interval: {name}/{key}")
        for suffix in ("_fit.png", "_eufspec.csv", "_flux_density.png", "_flux_density.csv"):
            require((directory / (name + suffix)).stat().st_size > 0,
                    f"Empty plot or data file: {name + suffix}")
        print(f"{name}: archived chi2/dof={fit['chi_square']:.5f}/{fit['dof']}, "
              f"L(1–10 keV)={lum['value']:.6e} erg/s")
    print("Saved fit/session/plot checks and luminosity arithmetic complete. XSPEC refits NOT VERIFIED.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("products", type=Path)
    directory = parser.parse_args().products
    events = {}
    for role in ("source", "background"):
        file = directory / f"{role}_events.fits"
        header, data = read_table(file, "EVENTS")
        require(int(header["OBS_ID"]) == 14946 and set(data["CCD_ID"]) == {7},
                f"Incorrect observation or CCD: {file.name}")
        _, gti = read_table(file, "GTI")
        start, stop = gti["START"], gti["STOP"]
        require(len(gti) > 0 and np.all(stop > start) and np.all(start[1:] >= stop[:-1]),
                f"Invalid GTI intervals: {file.name}")
        index = np.searchsorted(start, data["TIME"], side="right") - 1
        require(np.all(index >= 0) and np.all(data["TIME"] <= stop[index]),
                f"Events outside the GTI: {file.name}")
        exposure = float(np.sum(stop - start) * header["DTCOR"])
        require(math.isclose(exposure, header["EXPOSURE"], rel_tol=0, abs_tol=0.001),
                f"GTI/live exposure mismatch: {file.name}")
        events[role] = (exposure, data)
    ts, source = events["source"]
    tb, background = events["background"]
    hs, grouped = verify_pha(directory, "magnetar_grp.pi", source)
    hu, ungrouped = verify_pha(directory, "magnetar.pi", source)
    hb, back = verify_pha(directory, "magnetar_bkg.pi", background)
    require(np.array_equal(grouped["CHANNEL"], ungrouped["CHANNEL"]) and
            np.array_equal(grouped["CHANNEL"], back["CHANNEL"]), "PHA channel grids differ")
    for hdr, live in ((hs, ts), (hu, ts), (hb, tb)):
        require(math.isclose(hdr["EXPOSURE"], live, rel_tol=0, abs_tol=0.001),
                "PHA exposure differs from its event list")
    for hdr in (hs, hu):
        for key, file in (("BACKFILE", "magnetar_bkg.pi"), ("RESPFILE", "magnetar.rmf"),
                          ("ANCRFILE", "magnetar.corr.arf")):
            require(hdr[key] == file and (directory / file).is_file(),
                    f"Missing or incorrect association: {key}")
    _, ebounds = read_table(directory / "magnetar.rmf", "EBOUNDS")
    require(np.array_equal(ebounds["CHANNEL"], grouped["CHANNEL"]), "RMF channels differ from PHA")
    meta = json.loads((directory / "extraction_summary.json").read_text(encoding="utf-8"))
    alpha = (ts * hs["BACKSCAL"]) / (tb * hb["BACKSCAL"])
    areas = meta["areas_arcsec2"]
    require(math.isclose(alpha, areas["source"] / areas["background"], rel_tol=1e-7),
            "Exposure-scaled BACKSCAL ratio differs from the stored geometric area ratio")
    require(math.isclose(hu["BACKSCAL"], hs["BACKSCAL"], rel_tol=1e-7),
            "Grouped and ungrouped source BACKSCAL differ")
    ns, nb = (int(np.count_nonzero((ev["ENERGY"] >= 1000) & (ev["ENERGY"] < 10000)))
              for ev in (source, background))
    net = ns - alpha * nb
    sigma = math.sqrt(ns + alpha**2 * nb)
    print(f"GTI live times: source={ts:.9f} s, background={tb:.9f} s")
    print(f"All-channel counts: source={len(source)}, background={len(background)}")
    print(f"1–10 keV events: source={ns}, background={nb}, alpha={alpha:.9f}, "
          f"net={net:.5f} +/- {sigma:.5f}, rate={net/ts:.8f} +/- {sigma/ts:.8f} count/s")
    centres = (ebounds["E_MIN"] + ebounds["E_MAX"]) / 2
    mask = (centres >= 1) & (centres < 10)
    ch = ebounds["CHANNEL"][mask]
    print(f"1–10 keV PI-centre channels: {ch[0]}–{ch[-1]}, "
          f"counts={grouped['COUNTS'][mask].sum()}/{back['COUNTS'][mask].sum()}, "
          f"empty background channels={np.count_nonzero(back['COUNTS'][mask] == 0)}/{len(ch)}")
    grid = np.loadtxt(directory / "grouped_spectrum.csv", delimiter=",", skiprows=1)
    first, last = grid[0, 0] - grid[0, 1], grid[-1, 0] + grid[-1, 1]
    starts = np.flatnonzero(grouped["GROUPING"] != -1)
    ends = np.r_[starts[1:], len(grouped)]
    fitted = np.zeros(len(grouped), dtype=bool)
    for centre, half_width in grid[:, :2]:
        matches = np.flatnonzero(
            np.isclose(ebounds["E_MIN"][starts], centre - half_width, rtol=0, atol=5e-6)
            & np.isclose(ebounds["E_MAX"][ends - 1], centre + half_width, rtol=0, atol=5e-6)
        )
        require(len(matches) == 1, "Saved XSPEC bin is inconsistent with PHA grouping")
        start, stop = starts[matches[0]], ends[matches[0]]
        require(np.all(grouped["QUALITY"][start:stop] == 0), "Fit plot contains a bad group")
        require(grouped["COUNTS"][start:stop].sum() >= 50, "Fit plot contains a group under 50 counts")
        fitted[start:stop] = True
    print(f"Saved XSPEC plot: {len(grid)} groups spanning {first:.6f}–{last:.6f} keV")
    print(f"Nominal-band counts omitted by archived grouping/selection: "
          f"{grouped['COUNTS'][mask & ~fitted].sum()}/{back['COUNTS'][mask & ~fitted].sum()}")
    print("Full event selection, aspect histogram, DS9 and pileup run NOT VERIFIED.")
    verify_responses(directory)
    verify_fit_artifacts(directory)


if __name__ == "__main__":
    main()
