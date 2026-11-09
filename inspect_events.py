"""Measure aperture counts and plot archive or CIAO-reprocessed ACIS events."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from matplotlib.colors import PowerNorm
from matplotlib.patches import Circle
import warnings


ROOT = Path(__file__).resolve().parent
MAGNETAR = (266.4173708, -29.0082889)
SGRA = (266.416837, -29.007810)
MASKS = [
    (*SGRA, 1.5),
    (266.415746611, -29.008305206, 1.5),
    (266.416177055, -29.009646675, 1.5),
    (266.415826609, -29.006152406, 1.5),
    (266.416803127, -29.008626949, 1.0),
    (266.419632542, -29.007868751, 1.5),
]


def aperture_area(offset_masks):
    # Midpoint quadrature with 0.01 arcsec cells for the masked background.
    grid = np.arange(-7.995, 8, 0.01)
    xx, yy = np.meshgrid(grid, grid, sparse=True)
    radius = np.hypot(xx, yy)
    keep = (radius >= 2.5) & (radius < 8.0)
    for mx, my, size in offset_masks:
        keep &= np.hypot(xx - mx, yy - my) >= size
    return np.count_nonzero(keep) * 0.01**2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("events", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    with fits.open(args.events) as hdus:
        events = hdus["EVENTS"]
        header = events.header
        columns = [name.lower() for name in events.columns.names]
        xy_cols = [columns.index(axis) + 1 for axis in ("x", "y")]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(header, keysel=["pixel"], colsel=xy_cols, naxis=2)
        origin = wcs.all_world2pix([MAGNETAR], 1)[0]
        scale = np.abs(wcs.wcs.cdelt) * 3600
        x = (events.data["x"] - origin[0]) * scale[0]
        y = (events.data["y"] - origin[1]) * scale[1]
        energy = events.data["energy"] / 1000
        ccd = events.data["ccd_id"] == 7
        exposure = header["EXPOSURE"]
        offset_masks = [
            (*((wcs.all_world2pix([[ra, dec]], 1)[0] - origin) * scale), radius)
            for ra, dec, radius in MASKS
        ]
        sgra_x, sgra_y, _ = offset_masks[0]

        radius = np.hypot(x, y)
        source = ccd & (radius >= 0.75) & (radius < 2.0)
        background = ccd & (radius >= 2.5) & (radius < 8.0)
        for mx, my, size in offset_masks:
            background &= np.hypot(x - mx, y - my) >= size

        alpha = np.pi * (2.0**2 - 0.75**2) / aperture_area(offset_masks)
        band = (energy >= 1.0) & (energy < 10.0)
        ns = int(np.count_nonzero(source & band))
        nb = int(np.count_nonzero(background & band))
        summary = {
            "obsid": int(header["OBS_ID"]),
            "file": args.events.name,
            "caldb": header.get("CALDBVER"),
            "exposure_s": exposure,
            "frame_time_s": header["TIMEDEL"],
            "source_counts_1_10": ns,
            "background_counts_1_10": nb,
            "area_ratio": alpha,
            "net_counts_1_10": ns - alpha * nb,
            "net_counts_error_1_10": float(np.sqrt(ns + alpha**2 * nb)),
            "separation_from_reference_sgra_arcsec": float(np.hypot(sgra_x, sgra_y)),
        }

    image_band = ccd & (energy >= 1.0) & (energy < 8.0)
    edges = np.linspace(-10, 10, 81)
    image = np.histogram2d(x[image_band], y[image_band], bins=[edges, edges])[0]
    fig, ax = plt.subplots(figsize=(6.7, 5.1), layout="constrained")
    im = ax.imshow(image.T, origin="lower", extent=[-10, 10, -10, 10],
                   cmap="magma", norm=PowerNorm(0.38, vmin=0, vmax=75))
    for r in (0.75, 2):
        ax.add_patch(Circle((0, 0), r, fill=False, color="#65e4f0", lw=1.25))
    for r in (2.5, 8):
        ax.add_patch(Circle((0, 0), r, fill=False, color="#f5e276", lw=1.0, ls="--"))
    for mx, my, size in offset_masks:
        ax.add_patch(Circle((mx, my), size, fill=False, color="#eeeeee", lw=0.8, ls=":"))
    ax.annotate("SGR J1745−2900", xy=(-1.4, -1.4), xytext=(-9.3, -8.5), color="white",
                arrowprops={"arrowstyle": "-", "color": "white"}, fontsize=9)
    ax.annotate("Sgr A*", xy=(sgra_x, sgra_y), xytext=(3.8, 5.5), color="white",
                arrowprops={"arrowstyle": "-", "color": "white"}, fontsize=9)
    ax.set(xlabel="West offset from magnetar (arcsec)", ylabel="North offset (arcsec)",
           title="Chandra ACIS-S • ObsID 14946 • 1–8 keV")
    fig.colorbar(im, ax=ax, label="Counts per 0.25″ pixel", extend="max", shrink=0.82)
    fig.savefig(args.output / "field.png", dpi=240)
    plt.close(fig)

    bins = np.arange(1.0, 10.01, 0.25)
    s, _ = np.histogram(energy[source & band], bins)
    b, _ = np.histogram(energy[background & band], bins)
    net = s - alpha * b
    err = np.sqrt(s + alpha**2 * b)
    mid = (bins[1:] + bins[:-1]) / 2
    rate_scale = 1 / (exposure * np.diff(bins))
    fig, ax = plt.subplots(figsize=(6.9, 3.65), layout="constrained")
    ax.errorbar(mid, net * rate_scale, yerr=err * rate_scale,
                color="#14617d", fmt="o", ms=3.2, capsize=2,
                label="Background-subtracted aperture")
    ax.stairs(alpha * b * rate_scale, bins, color="#a66329", lw=1.35,
              label="Area-scaled background")
    ax.axhline(0, color="#8a8a8a", lw=0.7)
    ax.legend(frameon=False)
    ax.set(xlabel="Reconstructed event energy (keV)", ylabel="Count rate (counts s⁻¹ keV⁻¹)",
           xlim=(1, 10))
    fig.savefig(args.output / "counts.png", dpi=240)
    plt.close(fig)

    np.savetxt(args.output / "counts.csv", np.column_stack([mid, s, b, net, err, net * rate_scale,
                                                           err * rate_scale]),
               delimiter=",", comments="",
               header="energy_keV,source,background,net,net_poisson_sigma,rate_ct_s-1_keV-1,rate_sigma")
    (args.output / "event_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for filename, regions in (
        ("source.reg", [f'annulus({MAGNETAR[0]},{MAGNETAR[1]},0.75",2.0")']),
        ("background.reg", [f'annulus({MAGNETAR[0]},{MAGNETAR[1]},2.5",8.0")'] + [
            f'-circle({ra},{dec},{r}")' for ra, dec, r in MASKS
        ]),
    ):
        (args.output / filename).write_text(
            "# Region file format: DS9 version 4.1\nfk5\n" + "\n".join(regions) + "\n",
            encoding="ascii", newline="\n",
        )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
