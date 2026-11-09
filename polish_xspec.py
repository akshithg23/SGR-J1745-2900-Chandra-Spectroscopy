"""Supplementary checks for the saved SPA660 fits. Run in the HEASoft environment."""

import argparse
import csv
import json
import math
import os
from pathlib import Path

import numpy as np
from astropy.io import fits


DISTANCE_CM = 8.3 * 3.08567758149e21
LUMINOSITY_FACTOR = 4 * math.pi * DISTANCE_CM**2
FIT_NAME = "tbabs_bbody_free"


def restore_fit(products, row, arf=None, scale=1.0):
    from xspec import AllData, AllModels, Fit, Xset

    AllData.clear()
    AllModels.clear()
    Xset.restore(str(products / f"{FIT_NAME}.xcm"))
    Fit.query = "yes"
    Fit.nIterations = 200
    Fit.renorm("prefit")
    model = AllModels(1)
    model(1).values = row["nh_1e22_cm2"]["value"]
    model(2).values = row["shape"]["value"]
    model(3).values = row["norm"]["value"] / scale
    if arf is not None:
        AllData(1).response.arf = str(arf)
    return model


def check_statistic(row):
    from xspec import Fit

    if Fit.dof != row["dof"] or abs(Fit.statistic - row["chi_square"]) > 0.02:
        raise RuntimeError("The check changed the retained data or did not reproduce the saved fit")


def export_chart(products, output, row):
    from xspec import AllModels

    model = restore_fit(products, row)
    try:
        # ChaRT2 requires bin-integrated absorbed photon flux, without an ARF.
        for name, grid in (("chart_spectrum.dat", "1.0 10.0 900 lin"),
                           ("chart_4keV.dat", "3.99 4.01 1 lin")):
            AllModels.setEnergies(grid)
            edges = np.asarray(model.energies(1), dtype=float)
            photons = np.asarray(model.values(1), dtype=float)
            if len(edges) != len(photons) + 1 or np.any(photons <= 0):
                raise ValueError("Unusable ChaRT spectrum grid or photon flux")
            np.savetxt(output / name, np.column_stack((edges[:-1], edges[1:], photons)),
                       fmt="%.12e")
        AllModels.setEnergies("1.0 10.0 900 lin")
        total_photons = float(np.sum(model.values(1)))
        (output / "chart_total_photflux.txt").write_text(f"{total_photons:.12e}\n")
    finally:
        AllModels.setEnergies("reset")
    print("Wrote ChaRT2 spectra, using the saved absorbed blackbody without refitting.")
    print("Three columns: E_low [keV], E_high [keV], photons cm^-2 s^-1 per bin.")
    print("For a monochromatic 4-keV submission, use the total photon flux in")
    print(output / "chart_total_photflux.txt")


def sensitivity(products, output, row):
    from xspec import AllData, AllModels, Fit

    rows = []
    for scale in (1.0, 0.8, 1.2):
        target = output / f"magnetar_fpsf_{scale:.1f}.arf"
        with fits.open(products / "magnetar.corr.arf") as hdus:
            area = hdus["SPECRESP"].data
            area["SPECRESP"] *= scale
            area["PSF_FRAC"] *= scale
            if np.any((area["PSF_FRAC"] < 0) | (area["PSF_FRAC"] > 1)):
                raise ValueError("Scaled PSF fraction is outside [0,1]")
            hdus["SPECRESP"].header.add_history(f"Sensitivity test: f_PSF multiplied by {scale}")
            hdus.writeto(target, overwrite=True)
        model = restore_fit(products, row, target, scale)
        Fit.perform()
        check_statistic(row)
        nh, temperature, norm = (model(i).values[0] for i in (1, 2, 3))
        AllModels.calcFlux("1.0 10.0")
        absorbed = AllData(1).flux[0]
        model(1).values = 0.0
        try:
            AllModels.calcFlux("1.0 10.0")
            intrinsic = AllData(1).flux[0]
        finally:
            model(1).values = nh
        bolometric = norm * 1e39 * (8.3 / 10)**2
        radius_km = math.sqrt(bolometric / (4 * math.pi * 5.670374419e-5 *
                                          (temperature * 1.160451812e7)**4)) / 1e5
        rows.append(dict(fpsf_scale=scale, nh_1e22_cm2=nh, kT_keV=temperature,
                         norm=norm, chi_square=Fit.statistic, dof=Fit.dof,
                         absorbed_flux_1_10=absorbed, unabsorbed_flux_1_10=intrinsic,
                         luminosity_1_10=LUMINOSITY_FACTOR * intrinsic,
                         radius_km=radius_km))
    for item in rows:
        item["L_over_baseline"] = item["luminosity_1_10"] / rows[0]["luminosity_1_10"]
        item["R_over_baseline"] = item["radius_km"] / rows[0]["radius_km"]
        print(f"f scale={item['fpsf_scale']:.1f}, L={item['luminosity_1_10']:.6e} erg/s, "
              f"R={item['radius_km']:.5f} km, chi2/dof={item['chi_square']:.5f}/{item['dof']}")
    with (output / "psf_sensitivity.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def constrained_band(products, output, row):
    from xspec import AllData, AllModels, Fit, Model

    results = {}
    for kind in ("absorbed", "unabsorbed"):
        model = restore_fit(products, row)
        Fit.perform()
        check_statistic(row)
        nh, temperature, norm = (model(i).values[0] for i in (1, 2, 3))
        if kind == "unabsorbed":
            model(1).values = 0.0
        AllModels.calcFlux("2.26 6.69")
        seed = AllData(1).flux[0]
        AllModels.clear()
        if kind == "absorbed":
            model = Model("cflux*tbabs*bbody",
                          setPars=(2.26, 6.69, math.log10(seed), nh, temperature, norm))
            flux_index, fixed_indices = 3, (1, 2, 6)
        else:
            model = Model("tbabs*cflux*bbody",
                          setPars=(nh, 2.26, 6.69, math.log10(seed), temperature, norm))
            flux_index, fixed_indices = 4, (2, 3, 6)
        for index in fixed_indices:
            model(index).frozen = True
        Fit.perform()
        check_statistic(row)
        Fit.error(f"maximum 10.0 1.0 {flux_index}")
        check_statistic(row)
        low, high, status = model(flux_index).error
        centre = model(flux_index).values[0]
        valid = bool(status) and all(flag == "F" for flag in status) and low < centre < high
        result = dict(value=10**centre,
                      interval_68=[10**low, 10**high] if valid else None,
                      xspec_error_status=status, chi_square=Fit.statistic, dof=Fit.dof)
        results[kind + "_flux_2.26_6.69"] = result
        if kind == "unabsorbed":
            results["luminosity_2.26_6.69"] = dict(
                value=LUMINOSITY_FACTOR * result["value"],
                interval_68=[LUMINOSITY_FACTOR * x for x in result["interval_68"]]
                if valid else None)
    (output / "constrained_band.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
    print("Used the same retained groups. Only the cflux integration limits changed.")


def marx_fraction(paths, output, label):
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    import astropy.units as units

    target = SkyCoord(266.4173708, -29.0082889, unit="deg")
    total = annulus = 0
    for path in paths:
        with fits.open(path) as hdus:
            events = hdus["EVENTS"]
            xcol, ycol = (events.columns.names.index(name) + 1 for name in ("X", "Y"))
            wcs = WCS(events.header, hdus, keysel=["pixel"], colsel=[xcol, ycol])
            if wcs.pixel_n_dim != 2 or wcs.world_n_dim != 2:
                raise ValueError("The simulated EVENTS table did not provide a two-axis sky WCS")
            ra, dec = wcs.all_pix2world(events.data["X"], events.data["Y"], 1)
            radii = target.separation(SkyCoord(ra * units.deg, dec * units.deg)).arcsec
            # All detector-accepted simulated source events form the denominator.
            # No observed core counts, background events or finite radius are used.
            total += len(radii)
            annulus += int(np.count_nonzero((radii >= 0.75) & (radii < 2.0)))
    if not total:
        raise ValueError("No projected source events")
    fraction = annulus / total
    result = dict(label=label, total_simulated_source_events=total,
                  annular_events=annulus, annular_fraction=fraction,
                  binomial_standard_error=math.sqrt(fraction * (1 - fraction) / total),
                  input_files=[str(path) for path in paths])
    (output / f"marx_fraction_{label}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    print("Only a monochromatic 4-keV run is directly comparable with f_PSF(4 keV)=0.1077.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("chart", "sensitivity", "band", "fraction"))
    parser.add_argument("--products", type=Path, default=Path("products_verified"))
    parser.add_argument("--out", type=Path, default=Path("polish_checks"))
    parser.add_argument("--events", type=Path, nargs="+")
    parser.add_argument("--label", default="broadband")
    args = parser.parse_args()
    products, output = args.products.resolve(), args.out.resolve()
    if products == output or products in output.parents:
        raise ValueError("Use a supplementary output directory outside products_verified")
    output.mkdir(parents=True, exist_ok=True)
    if args.action == "fraction":
        if not args.events:
            parser.error("fraction requires --events")
        marx_fraction([path.resolve() for path in args.events], output, args.label)
        return
    from xspec import Xset

    row = json.loads((products / "fit_results.json").read_text())["fits"][FIT_NAME]
    os.chdir(products)
    Xset.openLog(str(output / f"{args.action}_xspec.log"))
    try:
        {"chart": export_chart, "sensitivity": sensitivity,
         "band": constrained_band}[args.action](products, output, row)
    finally:
        Xset.closeLog()


if __name__ == "__main__":
    main()
