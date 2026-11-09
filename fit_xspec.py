"""Fit an ObsID 14946 product directory with PyXspec."""

import argparse
import csv
import json
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from xspec import AllData, AllModels, Fit, Model, Plot, Spectrum, Xset


DISTANCE_KPC = 8.3
FROZEN_NH = 10.3  # XSPEC units of 10^22 cm^-2
KEV_TO_ERG = 1.602176634e-9
KEV_TO_K = 1.160451812e7
SIGMA_SB = 5.670374419e-5  # cgs
SPINDOWN = 4.9e33
LUMINOSITY_FACTOR = 4 * math.pi * (DISTANCE_KPC * 3.08567758149e21)**2


def check_responses():
    unc = fits.getdata("magnetar.arf", "SPECRESP")
    cor = fits.getdata("magnetar.corr.arf", "SPECRESP")
    if any(not np.allclose(unc[key], cor[key]) for key in ("ENERG_LO", "ENERG_HI")):
        raise ValueError("Corrected and uncorrected ARFs use different energy grids")
    energy = (cor["ENERG_LO"] + cor["ENERG_HI"]) / 2
    fraction = np.asarray(cor["PSF_FRAC"], dtype=float)
    band = (energy >= 1) & (energy <= 10) & (unc["SPECRESP"] > 0)
    frac = fraction[band]
    if not frac.size or np.any(~np.isfinite(frac) | (frac <= 0) | (frac > 1)):
        raise ValueError("Invalid 1–10 keV PSF fraction in the corrected ARF")
    if not np.allclose(cor["SPECRESP"][band], unc["SPECRESP"][band] * frac, rtol=1e-4, atol=0):
        raise ValueError("Corrected ARF does not agree with ARF times PSF_FRAC")
    if cor["ENERG_LO"][0] > 1 or cor["ENERG_HI"][-1] < 10:
        raise ValueError("ARF does not cover the complete 1–10 keV flux band")
    with fits.open("magnetar.rmf") as rmf:
        matrix = next(hdu.data for hdu in rmf[1:] if hdu.name.endswith("MATRIX"))
        for key in ("ENERG_LO", "ENERG_HI"):
            if matrix[key].shape != cor[key].shape or not np.allclose(matrix[key], cor[key]):
                raise ValueError("RMF and ARF use different incident-energy grids")
        if matrix["ENERG_LO"][0] > 1 or matrix["ENERG_HI"][-1] < 10:
            raise ValueError("RMF does not cover the complete 1–10 keV flux band")
    np.savetxt("psf_fraction.csv", np.column_stack([energy[band], frac, 1 / frac]), delimiter=",",
               header="energy_keV,annulus_psf_fraction,aperture_factor", comments="")
    f4 = float(np.interp(4, energy[band], frac))
    return {
        "fraction_at_4_keV": f4,
        "factor_at_4_keV": 1 / f4,
        "min_fraction_1_10": float(frac.min()),
        "max_fraction_1_10": float(frac.max()),
    }


def interval(parameter):
    low, high, status = parameter.error
    value = parameter.values[0]
    valid = status == "FFFFFFFFF" and low < value < high
    return {
        "value": value,
        "interval_68": [low, high] if valid else None,
        "xspec_error_status": status,
    }


def save_session(stem):
    path = Path(f"{stem}.xcm")
    path.unlink(missing_ok=True)
    Xset.save(str(path))


def plot_spectrum(stem, fitted=False):
    Plot("data", "delchi") if fitted else Plot("data")
    if fitted:
        fig, (top, bottom) = plt.subplots(
            2, 1, figsize=(7, 5.4), sharex=True, layout="constrained",
            gridspec_kw={"height_ratios": [3, 1]},
        )
    else:
        fig, top = plt.subplots(figsize=(7, 3.8), layout="constrained")
        bottom = top

    energy, counts = Plot.x(1, 1), Plot.y(1, 1)
    width, uncertainty = Plot.xErr(1, 1), Plot.yErr(1, 1)
    top.errorbar(energy, counts, xerr=width, yerr=uncertainty, fmt="o", ms=3,
                 color="#14617d", label="Chandra data")
    top.set_ylabel("Counts s⁻¹ keV⁻¹")
    if fitted:
        top.plot(energy, Plot.model(1, 1), color="#a66329", label="Folded model")
        top.legend(frameon=False)
        bottom.errorbar(Plot.x(1, 2), Plot.y(1, 2), xerr=Plot.xErr(1, 2),
                        yerr=Plot.yErr(1, 2), fmt="o", ms=3, color="#14617d")
        bottom.axhline(0, color="0.5", lw=1)
        bottom.set_ylabel("Residual (σ)")
    bottom.set(xlabel="Energy (keV)", xlim=(1, 10))
    fig.savefig(f"{stem}.png", dpi=240)
    plt.close(fig)
    if not fitted:
        np.savetxt(f"{stem}.csv", np.column_stack([energy, width, counts, uncertainty]),
                   delimiter=",", comments="",
                   header="energy_keV,half_bin_width_keV,counts_per_s_keV,sigma_counts_per_s_keV")


def save_flux_density(model, nh, stem):
    # eufspec points depend on the fitted model. Convert keV-based energy
    # flux density to erg-based units, keeping the plot bins XSPEC used.
    Plot("eufspec")
    unfolded = np.column_stack([Plot.x(), Plot.xErr(), Plot.y(), Plot.yErr(), Plot.model()])
    unfolded[:, 2:] *= KEV_TO_ERG
    np.savetxt(f"{stem}_eufspec.csv", unfolded, delimiter=",", comments="",
               header="energy_keV,half_bin_width_keV,unfolded_FE_cgs_per_keV,"
                      "unfolded_sigma_cgs_per_keV,absorbed_model_FE_cgs_per_keV")

    # Evaluate the fitted continuum on a fine energy grid independently of
    # grouping. Model.values returns photons integrated in each model bin.
    AllModels.setEnergies("1.0 10.0 1800 lin")
    try:
        edges = np.asarray(model.energies(1))
        energy = (edges[:-1] + edges[1:]) / 2
        factor = energy * KEV_TO_ERG / np.diff(edges)
        absorbed = np.asarray(model.values(1)) * factor
        model(1).values = 0.0
        intrinsic = np.asarray(model.values(1)) * factor
    finally:
        model(1).values = nh
        AllModels.setEnergies("reset")

    np.savetxt(f"{stem}_flux_density.csv", np.column_stack([energy, absorbed, intrinsic]),
               delimiter=",", comments="",
               header="energy_keV,absorbed_erg_cm-2_s-1_keV-1,unabsorbed_erg_cm-2_s-1_keV-1")
    fig, ax = plt.subplots(figsize=(7, 3.8), layout="constrained")
    ax.errorbar(unfolded[:, 0], unfolded[:, 2], xerr=unfolded[:, 1], yerr=unfolded[:, 3],
                fmt="o", ms=2.5, color="0.5", alpha=0.65, label="XSPEC unfolded data")
    ax.plot(energy, absorbed, color="#14617d", label="Absorbed model")
    ax.plot(energy, intrinsic, color="#a66329", label="Intrinsic model")
    ax.set(xlabel="Photon energy (keV)", ylabel="Energy flux density (erg cm⁻² s⁻¹ keV⁻¹)",
           xlim=(1, 10))
    ax.legend(frameon=False)
    fig.savefig(f"{stem}_flux_density.png", dpi=240)
    plt.close(fig)

    def band_peak(density):
        peak = int(np.argmax(density))
        if peak in (0, len(energy) - 1):
            return {"energy_keV": float(edges[0 if peak == 0 else -1]), "at_band_edge": True}
        return {"energy_keV": float(energy[peak]), "at_band_edge": False}

    return {
        "absorbed_at_4_keV": float(np.interp(4.0, energy, absorbed)),
        "unabsorbed_at_4_keV": float(np.interp(4.0, energy, intrinsic)),
        "unit": "erg cm^-2 s^-1 keV^-1",
        "absorbed_band_peak": band_peak(absorbed),
        "unabsorbed_band_peak": band_peak(intrinsic),
    }


def cflux_fit(absorption, continuum, values, fixed, unabsorbed, seed, baseline_chi):
    nh, shape, norm = values
    if seed <= 0 or norm <= 0:
        raise RuntimeError("Cannot initialize cflux with a nonpositive flux or continuum normalization")
    log_flux = math.log10(seed)
    AllModels.clear()
    if unabsorbed:
        model = Model(f"{absorption}*cflux*{continuum}", setPars=(nh, 1, 10, log_flux, shape, norm))
        nh_index, emin, emax, flux_index = 1, 2, 3, 4
    else:
        model = Model(f"cflux*{absorption}*{continuum}", setPars=(1, 10, log_flux, nh, shape, norm))
        emin, emax, flux_index, nh_index = 1, 2, 3, 4
    for index in (emin, emax, 6):
        model(index).frozen = True
    model(nh_index).frozen = fixed
    Fit.perform()
    Fit.error(f"maximum 10.0 1.0 {flux_index}")
    if not math.isfinite(Fit.statistic):
        raise RuntimeError("cflux returned a nonfinite fit statistic")
    estimate = interval(model(flux_index))
    estimate["value"] = 10**estimate["value"]
    if estimate["interval_68"]:
        estimate["interval_68"] = [10**bound for bound in estimate["interval_68"]]
    estimate["chi_square"] = Fit.statistic
    estimate["delta_chi_from_baseline"] = Fit.statistic - baseline_chi
    if abs(estimate["delta_chi_from_baseline"]) > 0.05:
        raise RuntimeError("The cflux fit differs from the baseline χ² by more than 0.05")
    return estimate


def run_fit(spectrum, absorption, continuum, initial, fixed):
    stem = f"{absorption}_{continuum}_{'fixed' if fixed else 'free'}"
    AllModels.clear()
    model = Model(f"{absorption}*{continuum}", setPars=initial)
    if fixed:
        model(1).values = FROZEN_NH
        model(1).frozen = True
    Fit.perform()
    if Fit.dof <= 0 or not math.isfinite(Fit.statistic):
        raise RuntimeError(f"Unusable baseline fit: {stem}")
    Fit.error("maximum 10.0 1.0 2-3" if fixed else "maximum 10.0 1.0 1-3")
    values = tuple(model(i).values[0] for i in range(1, 4))
    nh, shape, norm = values
    chi_square = Fit.statistic
    result = {
        "model": model.expression,
        "nh_1e22_cm2": {"value": nh, "interval_68": None, "fixed": True} if fixed else interval(model(1)),
        "shape": interval(model(2)),
        "norm": interval(model(3)),
        "chi_square": chi_square,
        "dof": Fit.dof,
        "reduced_chi_square": chi_square / Fit.dof,
        "chi_square_null_probability": Fit.nullhyp,
    }
    save_session(stem)
    plot_spectrum(f"{stem}_fit", fitted=True)
    result["flux_density"] = save_flux_density(model, nh, stem)

    # These ordinary fluxes initialize the corresponding cflux parameters.
    AllModels.calcFlux("1.0 10.0")
    seed_abs = spectrum.flux[0]
    try:
        model(1).values = 0.0
        AllModels.calcFlux("1.0 10.0")
        seed_unabs = spectrum.flux[0]
    finally:
        model(1).values = nh
    result["absorbed_flux_1_10"] = cflux_fit(
        absorption, continuum, values, fixed, False, seed_abs, chi_square)
    intrinsic = cflux_fit(absorption, continuum, values, fixed, True, seed_unabs, chi_square)
    result["unabsorbed_flux_1_10"] = intrinsic
    result["luminosity_1_10"] = {
        "value": LUMINOSITY_FACTOR * intrinsic["value"],
        "interval_68": ([LUMINOSITY_FACTOR * flux for flux in intrinsic["interval_68"]]
                        if intrinsic["interval_68"] else None),
    }
    result["luminosity_over_spindown"] = result["luminosity_1_10"]["value"] / SPINDOWN
    if continuum == "bbody":
        if norm <= 0 or shape <= 0:
            raise ValueError("Nonphysical blackbody temperature or normalization")
        bolometric_luminosity = 1e39 * norm * (DISTANCE_KPC / 10)**2
        radius = math.sqrt(bolometric_luminosity / (4 * math.pi * SIGMA_SB * (shape * KEV_TO_K)**4))
        result["blackbody_radius_km"] = radius / 1e5
        result["intrinsic_FE_peak_keV"] = 2.821439 * shape
        result["intrinsic_EFE_peak_keV"] = 3.920690 * shape
    return stem, result


def write_comparison(results):
    reference_flux = 1.6e35 / LUMINOSITY_FACTOR
    with Path("rea_comparison.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["fit", "quantity", "unit", "measured", "lower_68", "upper_68",
                         "Rea_et_al_or_derived", "percent_difference", "reference_scope"])
        for name, row in results.items():
            blackbody = "_bbody_" in name
            entries = [
                ("N_H", "1e22 cm^-2", row["nh_1e22_cm2"], 10.3 if blackbody else 17.3,
                 "three ACIS-S observations"),
                ("kT" if blackbody else "Gamma", "keV" if blackbody else "dimensionless",
                 row["shape"], 0.88 if blackbody else 4.2, "ObsID 14946"),
                ("absorbed flux", "erg cm^-2 s^-1", row["absorbed_flux_1_10"], 6.3e-12,
                 "ObsID 14946 Table 1"),
                ("unabsorbed flux", "erg cm^-2 s^-1", row["unabsorbed_flux_1_10"], reference_flux,
                 "derived from published blackbody luminosity" if blackbody else "blackbody benchmark only"),
                ("luminosity", "erg s^-1", row["luminosity_1_10"], 1.6e35,
                 "ObsID 14946 blackbody at 8.3 kpc" if blackbody else "blackbody benchmark only"),
            ]
            for quantity, unit, measured, ref, scope in entries:
                lo, hi = measured["interval_68"] or ("", "")
                value = measured["value"]
                pct = "" if measured.get("fixed") else 100 * (value - ref) / ref
                if measured.get("fixed"):
                    scope += "; fitted column was fixed by design"
                writer.writerow([name, quantity, unit, value, lo, hi, ref, pct, scope])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--products", type=Path, default=Path(__file__).resolve().parent / "products")
    parser.add_argument("--overwrite-results", action="store_true",
                        help="Allow replacement of existing fit products in this run directory")
    parser.add_argument("--rea-crosscheck", action="store_true",
                        help="also fit free-column phabs models with angr/bcmc for comparison to Rea et al.")
    args = parser.parse_args()
    if (args.products / "fit_results.json").exists() and not args.overwrite_results:
        parser.error("Fit outputs already exist. Use a separate run directory or --overwrite-results")
    os.chdir(args.products)
    Path("fit_results.json").unlink(missing_ok=True)
    response = check_responses()
    Xset.openLog("xspec.log")
    try:
        Xset.abund = "angr"
        Xset.xsect = "bcmc"  # Honored by phabs. TBabs supplies its own cross-sections.
        Xset.addModelString("TBABSVERSION", "2")
        spectrum = Spectrum("magnetar_grp.pi", backFile="magnetar_bkg.pi",
                            respFile="magnetar.rmf", arfFile="magnetar.corr.arf")
        Plot.device = "/null"
        Plot.xAxis = "keV"
        AllData.ignore("bad")
        spectrum.ignore("**-1.0 10.0-**")
        Fit.statMethod = "chi"
        Fit.weight = "standard"
        # Keep the fitted normalization when evaluating the unabsorbed model.
        Fit.renorm("prefit")
        Fit.query = "yes"
        Fit.nIterations = 200
        plot_spectrum("grouped_spectrum")
        results = {}
        absorption_models = ("tbabs", "phabs") if args.rea_crosscheck else ("tbabs",)
        for absorption in absorption_models:
            for continuum, initial in (("bbody", (10.3, 0.9, 2.0e-4)),
                                       ("powerlaw", (17.3, 4.2, 0.1))):
                stem, result = run_fit(spectrum, absorption, continuum, initial, fixed=False)
                results[stem] = result
                if absorption == "tbabs":
                    start = (FROZEN_NH, result["shape"]["value"], result["norm"]["value"])
                    stem, result = run_fit(spectrum, absorption, continuum, start, fixed=True)
                    results[stem] = result

        comparison = {}
        for label in ("free", "fixed"):
            bb, pl = (results[f"tbabs_{component}_{label}"] for component in ("bbody", "powerlaw"))
            comparison[label] = bb["chi_square"] - pl["chi_square"]
        write_comparison(results)
        Path("fit_results.json").write_text(json.dumps({
            "xspec_version": Xset.version,
            "abund": "angr",
            "xsect": "bcmc (phabs only)",
            "tbabs_version": 2,
            "energy_band_keV": [1, 10],
            "distance_kpc": DISTANCE_KPC,
            "aperture": response,
            "fits": results,
            "chi_square_blackbody_minus_powerlaw": comparison,
        }, indent=2, allow_nan=False) + "\n")
    finally:
        Xset.closeLog()


if __name__ == "__main__":
    main()
