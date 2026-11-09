"""Build the SPA660 report from the completed SciServer results archive."""

import argparse
import csv
from html import escape
from io import BytesIO, StringIO
import json
import math
from pathlib import Path
from zipfile import ZipFile

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Image, PageBreak, SimpleDocTemplate, Spacer, Table

from build_report import footer, href, p, styles, table as base_table

# Keep the complete bibliography legible while allowing it to fit after the
# conclusion without creating a references-only page.
styles["ref"] = ParagraphStyle(
    "FinalReference", fontName="Report", fontSize=9.5, leading=11.5,
    textColor="#202833", spaceAfter=2,
)


FIT_LABELS = {
    "tbabs_bbody_free": "TBabs blackbody, free",
    "tbabs_bbody_fixed": "TBabs blackbody, fixed",
    "tbabs_powerlaw_free": "TBabs power law, free",
    "tbabs_powerlaw_fixed": "TBabs power law, fixed",
    "phabs_bbody_free": "phabs blackbody, free",
    "phabs_powerlaw_free": "phabs power law, free",
}


def table(rows, widths):
    result = base_table(rows, widths)
    result.setStyle([("TOPPADDING", (0, 0), (-1, -1), 4.5),
                     ("BOTTOMPADDING", (0, 0), (-1, -1), 4.5)])
    return result


def interval(entry, scale=1.0, digits=2):
    value = entry["value"] / scale
    bounds = entry["interval_68"]
    if entry.get("fixed"):
        return f"{value:.1f} (fixed)"
    if bounds is None:
        return f"{value:.{digits}f} (interval unresolved)"
    low, high = (bound / scale for bound in bounds)
    return f"{value:.{digits}f} [{low:.{digits}f}, {high:.{digits}f}]"


def graphic(data, width):
    item = Image(BytesIO(data))
    item.drawWidth = width
    item.drawHeight = item.imageHeight * width / item.imageWidth
    return item


def comparison(label, entry, reference, ref_text, scale=1.0):
    delta = 100 * (entry["value"] / reference - 1)
    return [label, interval(entry, scale), ref_text, f"{delta:+.1f}".replace("-", "−")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zip_path", type=Path)
    parser.add_argument("pdf_path", type=Path)
    parser.add_argument("--code-url", help="Public-facing URL for the project code")
    args = parser.parse_args()
    zip_path, pdf_path = args.zip_path, args.pdf_path
    with ZipFile(zip_path) as archive:
        def read(name):
            return archive.read(f"products_verified/{name}")

        summary = json.loads(read("event_check/event_summary.json"))
        extraction = json.loads(read("extraction_summary.json"))
        results = json.loads(read("fit_results.json"))
        fits = results["fits"]
        bb = fits["tbabs_bbody_free"]
        pl = fits["tbabs_powerlaw_free"]
        pl_fixed = fits["tbabs_powerlaw_fixed"]
        phabs_bb = fits["phabs_bbody_free"]
        phabs_pl = fits["phabs_powerlaw_free"]
        aperture = results["aperture"]
        field = graphic(read("event_check/field.png"), 480)
        count_plot = graphic(read("event_check/counts.png"), 485)
        fit_pl = graphic(read("tbabs_powerlaw_free_fit.png"), 240)
        fit_bb = graphic(read("tbabs_bbody_free_fit.png"), 240)
        energy_plot = graphic(read("tbabs_bbody_free_flux_density.png"), 447)
        density = list(csv.reader(StringIO(read("tbabs_bbody_free_flux_density.csv").decode())))[1:]
        group_rows = list(csv.reader(StringIO(read("grouped_spectrum.csv").decode())))[1:]

    lo = float(group_rows[0][0]) - float(group_rows[0][1])
    hi = float(group_rows[-1][0]) + float(group_rows[-1][1])
    dist = results["distance_kpc"]
    ref_flux = 1.6e35 / (4 * math.pi * (dist * 3.08567758149e21)**2)
    live = extraction["exposure_s"]
    net = summary["net_counts_1_10"]
    sigma = summary["net_counts_error_1_10"]
    alpha = summary["area_ratio"]
    net_text = f"{net:,.0f} ± {sigma:.0f}"
    dx = results["chi_square_blackbody_minus_powerlaw"]["free"]
    reference_link = "https://doi.org/10.1088/2041-8205/775/2/L34"
    models = [
        [name, interval(fits[key]["nh_1e22_cm2"], digits=1),
         interval(fits[key]["shape"]), f"{fits[key]['chi_square']:.2f} / {fits[key]['dof']}"]
        for key, name in FIT_LABELS.items()
    ]
    fluxes = [
        [FIT_LABELS[key], interval(fits[key]["absorbed_flux_1_10"], 1e-11),
         interval(fits[key]["unabsorbed_flux_1_10"], 1e-11),
         interval(fits[key]["luminosity_1_10"], 1e35)]
        for key in FIT_LABELS if key.startswith("tbabs")
    ]
    comparison_rows = [
        comparison("Blackbody N<sub>H</sub> (10²² cm⁻²)", phabs_bb["nh_1e22_cm2"], 10.3, "10.3 ± 0.4"),
        comparison("Power-law N<sub>H</sub> (10²² cm⁻²)", phabs_pl["nh_1e22_cm2"], 17.3, "17.3 ± 0.4"),
        comparison("Blackbody kT (keV)", phabs_bb["shape"], 0.88, "0.88 ± 0.02"),
        comparison("Power-law Γ", phabs_pl["shape"], 4.2, "4.2 ± 0.1"),
        comparison("Absorbed flux (10⁻¹² erg cm⁻² s⁻¹)", phabs_bb["absorbed_flux_1_10"], 6.3e-12, "6.3 ± 0.2", 1e-12),
        comparison("Intrinsic flux (10⁻¹¹ erg cm⁻² s⁻¹)", phabs_bb["unabsorbed_flux_1_10"], ref_flux, "≈ 1.9 (derived)", 1e-11),
        comparison("Luminosity (10³⁵ erg s⁻¹)", phabs_bb["luminosity_1_10"], 1.6e35, "1.6", 1e35),
    ]
    density_rows = []
    grid = [float(row[0]) for row in density]
    for energy in (2, 4, 6, 8):
        index = next(i for i, value in enumerate(grid) if value > energy)
        low, high = density[index - 1], density[index]
        weight = (energy - float(low[0])) / (float(high[0]) - float(low[0]))
        values = [float(low[j]) * (1 - weight) + float(high[j]) * weight for j in (1, 2)]
        density_rows.append([str(energy), *(f"{value/1e-12:.3g}" for value in values)])

    panels = Table(
        [[p("<b>Blackbody with free column</b>", "caption"), p("<b>Power law with free column</b>", "caption")],
         [fit_bb, fit_pl]],
        colWidths=[246, 246],
    )
    panels.setStyle([("LEFTPADDING", (0, 0), (-1, -1), 3),
                     ("RIGHTPADDING", (0, 0), (-1, -1), 3)])

    story = [
        p("SPA660 ASTROPHYSICS PROJECT  /  PROJECT 8", "kicker"),
        p("X-ray spectroscopy of the<br/>Galactic Centre magnetar", "title"),
        p("SGR J1745−2900 in Chandra ObsID 14946", "subtitle"),
        p("SPA660 Project 8 group"),
        Spacer(1, 6),
        p("Abstract", "h1"),
        p(f"We reprocessed an 18.203 ks Chandra ACIS-S observation of SGR J1745−2900. "
          f"The 0.75–2.0 arcsec annulus yielded {net_text} net counts in 1–10 keV. "
          "We applied an energy-dependent annular PSF correction. Free-column blackbody and "
          "power-law fits gave χ²/dof "
          f"of {bb['chi_square']:.2f}/24 and {pl['chi_square']:.2f}/24. "
          "The blackbody remained adequate with the column fixed at the published value. Its "
          "free-column fit gave kT = 0.81 keV, R<sub>BB</sub> = 2.99 km, and "
          "L<sub>X</sub> = 4.70 × 10³⁵ erg s⁻¹ (2.9 times Rea et al.; radius 2.1 times theirs). "
          "The matched phabs N<sub>H</sub> ≈ 1.1 × 10²³ cm⁻² was within 7% of Rea et al.; absorbed flux "
          "was 2.4 times Rea et al. and luminosity was 96 times the spin-down power. The 1–10 keV "
          "flux and luminosity remained extrapolation-dependent because groups covered only "
          "2.26–6.69 keV."),
        p("1  Target and observation", "h1"),
        p("SGR J1745−2900 is a transient magnetar near the Galactic Centre. Its published "
          "position lies 2.4 ± 0.3 arcsec from Sgr A* [1]. Chandra resolves the two positions "
          "in a field with strong diffuse X-ray emission."),
        table([
            ["<b>Observation property</b>", "<b>Value</b>"],
            ["Chandra observation", "ObsID 14946 · 2 July 2013"],
            ["Detector and mode", "ACIS-S3 (CCD 7) · FAINT · 1/8 subarray"],
            ["CCD 7 exposure and frame time", f"{live/1000:.3f} ks · {summary['frame_time_s']:.5f} s"],
            ["Reprocessing software and calibration", "CIAO 4.17.0 · CALDB 4.12.4 (processing history)"],
        ], [199, 286]),
        Spacer(1, 8),
        p("2  Reprocessing and event selection", "h1"),
        p("We reprocessed the archive products [2] with <b>chandra_repro</b>. The CCD 7 "
          "aspect-histogram exposure agreed with the GTI × live-time exposure to within 0.001 s."),
        p("We selected events in sky coordinates around RA 266.4173708° and Dec −29.0082889° "
          "(J2000); the adopted position was confirmed in DS9. The 0.75–2.0 arcsec source annulus "
          "excludes the bright core to reduce pileup [1]. The local 2.5–8.0 arcsec background "
          "annulus used six masks for Sgr A* and nearby catalogue sources [3]."),
        p("<b>Pileup.</b> <b>pileup_map</b> maxima were 0.154 counts/frame in the excluded core "
          "and 0.0627 counts/frame in the source annulus. The 0.1 counts/frame value "
          "is a rough guide to about 5% pileup [7]. The core exceeded this guide while the annulus "
          "remained below it, supporting the 0.75 arcsec inner exclusion; no pileup fraction was measured."),
        PageBreak(),
        p("3  Galactic Centre image", "h1"),
        field,
        p("<b>Figure 1.</b> Reprocessed Chandra events from 1–8 keV. Solid cyan lines mark the "
          "source annulus, dashed yellow lines bound the background annulus, and dotted white "
          "circles mark exclusions from that background. The colour stretch saturates the core "
          "for visibility. The axes give offsets from the adopted magnetar position.", "caption"),
        p(f"The bright core lies at the published magnetar position. The adopted coordinates "
          f"place Sgr A* {summary['separation_from_reference_sgra_arcsec']:.2f} arcsec to the northwest [1,3]."),
        p("Aperture photometry", "h2"),
        table([
            ["<b>Quantity</b>", "<b>1–10 keV result</b>"],
            ["Source annulus events", f"{summary['source_counts_1_10']:,}"],
            ["Masked background annulus events", f"{summary['background_counts_1_10']:,}"],
            ["Source-to-background area ratio", f"{alpha:.5f}"],
            ["Net source counts", net_text],
            ["Net annular count rate", f"{net/live:.4f} ± {sigma/live:.4f} counts s⁻¹"],
        ], [275, 210]),
        p("The net counts are N<sub>s</sub> − αN<sub>b</sub>, where α is the geometric area ratio. "
          "The counting uncertainty is √(N<sub>s</sub> + α²N<sub>b</sub>). The rate describes the "
          "observed annulus before the point spread function correction.", "caption"),
        PageBreak(),
        p("4  Spectrum extraction and calibration", "h1"),
        count_plot,
        p("<b>Figure 2.</b> Background-subtracted detector count-rate density in 0.25 keV bins. "
          "The orange histogram is the area-scaled background. The peak detector bin is "
          "3.50–3.75 keV and the counts decline above 6 keV. This is a detector-space reference plot.", "caption"),
        p("PI spectrum and response", "h2"),
        p("The region-filtered extraction returned an empty spectrum, so we selected the source "
          "and background events from the reprocessed table in physical coordinates and ran CIAO "
          "<b>dmextract</b> on the resulting event lists. We verified the PI counts channel by channel "
          "against the selected events. Both spectra use the checked CCD 7 exposure. "
          "Their BACKSCAL ratio is the geometric area ratio after applying the background masks. "
          "<b>specextract</b> generated the unweighted ARF and RMF using the explicit aspect histogram [4]."),
        p("Point spread function correction", "h2"),
        p("We used CIAO’s CALDB REEF calibration to obtain the radial encircled count fractions [5]. "
          "For each ARF energy bin we subtracted the fraction inside 0.75 arcsec from the fraction "
          "inside 2.0 arcsec. The corrected ARF multiplies the effective area by this annular fraction. "
          f"At 4 keV the fraction is {aperture['fraction_at_4_keV']:.4f}, giving a reciprocal "
          f"correction factor of {aperture['factor_at_4_keV']:.2f}. Across 1–10 keV it ranges from "
          f"{aperture['min_fraction_1_10']:.4f} to {aperture['max_fraction_1_10']:.4f}. "
          "An independent CIAO arfcorr check gave 0.1092 at 4 keV and 0.07880–0.1258 over 1–10 keV."),
        p("Grouping and usable energy range", "h2"),
        p(f"CIAO <b>dmgroup</b> grouped the source PI spectrum to at least 50 counts per good group. "
          f"XSPEC selected 1–10 keV and rejected bad-quality bins. The {len(group_rows)} retained "
          f"groups actually cover {lo:.2f}–{hi:.2f} keV and contain 50–64 source counts per group. "
          "The broad low-energy group crosses the lower boundary, and the final high-energy group "
          "falls below the count threshold. The quoted 1–10 keV flux therefore depends on the "
          "model outside the retained groups."),
        PageBreak(),
        p("5  XSPEC fits and reference comparison", "h1"),
        p(f"We used XSPEC {results['xspec_version'][1]} with χ² and standard weighting [6]. "
          "The TBabs fits use the version 2 absorption prescription and the XSPEC session abundance "
          "setting <b>angr</b>. The two <b>phabs</b> cross-checks use <b>angr</b> and <b>bcmc</b> "
          "to match the absorption settings in Rea et al. [1]. Brackets give one-parameter "
          "68% profile intervals from Δχ² = 1."),
        table([
            ["<b>Continuum and column</b>", "<b>N<sub>H</sub> (10²² cm⁻²)</b>", "<b>kT (keV) or Γ</b>", "<b>χ²/dof</b>"],
            *models,
        ], [154, 125, 129, 77]),
        p("<b>Table 1.</b> Both fixed-N<sub>H</sub> fits used N<sub>H</sub> = 10.3 × 10²² cm⁻², "
          "the published blackbody column.", "caption"),
        table([
            ["<b>TBabs fit</b>", "<b>Absorbed flux</b>", "<b>Intrinsic flux</b>", "<b>L<sub>X</sub></b>"],
            *fluxes,
        ], [147, 105, 122, 111]),
        p("<b>Table 2.</b> Fluxes are in 10⁻¹¹ erg cm⁻² s⁻¹ and luminosities in 10³⁵ erg s⁻¹ "
          "over 1–10 keV. We used <b>cflux</b> on the full model and on the continuum alone, "
          "then calculated L<sub>X</sub> = 4πd²F<sub>unabs</sub> at 8.3 kpc. "
          "The fixed power law’s intervals are conditional on a poor fit.", "caption"),
        table([
            ["<b>Quantity and unit</b>", "<b>Our phabs fit</b>", "<b>Rea et al.</b>", "<b>Δ (%)</b>"],
            *comparison_rows,
        ], [188, 120, 105, 72]),
        p("<b>Table 3.</b> Δ = 100 × (fit − reference)/reference. Published columns are from "
          "three ACIS-S observations. Temperatures and indices are for ObsID 14946. "
          "The flux and luminosity comparison uses the free blackbody. The reference intrinsic "
          "flux is derived from the reported luminosity at 8.3 kpc [1].", "caption"),
        PageBreak(),
        p("6  Spectral shape and residuals", "h1"),
        panels,
        p("<b>Figure 3.</b> Free-column TBabs fits. The upper panels show the background-subtracted "
          "data and response-folded model. The lower panels show residuals in σ. Both models show "
          "a positive residual near 3.4 keV and a negative residual near 4 keV.", "caption"),
        energy_plot,
        p("<b>Figure 4.</b> Energy flux density for the free-column blackbody. The plotted "
          "unfolded data depend on the fitted model. The curves show the absorbed and intrinsic model.", "caption"),
        table([
            ["<b>Energy (keV)</b>", "<b>Absorbed F<sub>E</sub></b>", "<b>Intrinsic F<sub>E</sub></b>"],
            *density_rows,
        ], [157, 164, 164]),
        p("<b>Table 4.</b> Model energy flux densities in 10⁻¹² erg cm⁻² s⁻¹ keV⁻¹. "
          "The full grid and XSPEC unfolded values are provided in the results archive.", "caption"),
        PageBreak(),
        p("7  Physical interpretation", "h1"),
        p(f"<b>Continuum and source nature.</b> Both free-column fits were adequate (Δχ² = {dx:.2f}), "
          f"leaving the continuum model-dependent; at the published column the power law fit poorly "
          f"(χ²/dof = {pl_fixed['chi_square']:.2f}/25), overpredicting both band ends. Blackbody "
          f"K = {bb['norm']['value'] * 1e4:.2f} × 10⁻⁴ gave R<sub>BB</sub> = "
          f"{bb['blackbody_radius_km']:.2f} km from L<sub>bol</sub> = 4πR²σT⁴ at 0.81 keV and "
          "8.3 kpc (2.1 times Rea et al.’s ≈1.4 km). With kT = 0.81 versus 0.88 keV, R²T⁴ was "
          "about 3.3, near the 2.9-fold luminosity ratio; the radius was consistent with a hot region."),
        p("<b>Absorption and energy scale.</b> Free blackbody N<sub>H</sub> = 1.12 × 10²³ cm⁻² "
          "indicated heavy Galactic Centre absorption. At 4 keV, absorbed/intrinsic F<sub>E</sub> was "
          "4.63 × 10⁻¹²/1.00 × 10⁻¹¹ erg cm⁻² s⁻¹ keV⁻¹, peaking at 4.04/2.29 keV. The power law "
          "rose to the lower energy edge without an interior peak."),
        p("<b>Luminosity.</b> Blackbody L<sub>X</sub> = 4.70 × 10³⁵ erg s⁻¹ "
          "(68% interval 4.09–5.46 × 10³⁵) was 96 times the 4.9 × 10³³ erg s⁻¹ spin-down power [1]; "
          "rotation alone could not sustain it. The steep power law gave about 10³⁷ erg s⁻¹ by "
          "extrapolation into the absorbed soft band. Thus L<sub>X</sub> over 1–10 keV depended on "
          "extrapolation beyond groups at 2.26–6.69 keV."),
        p("<b>Agreement and limitations.</b> Phabs blackbody N<sub>H</sub> = 1.10 × 10²³ cm⁻² "
          "was within 7% of Rea; kT = 0.81 keV (8% lower), L<sub>X</sub> = 4.70 × 10³⁵ erg s⁻¹ "
          "(2.9×), and R<sub>BB</sub> = 2.99 km (2.1×). TBabs blackbody/power-law χ²/dof were "
          "23.06/24 and 20.88/24; phabs L<sub>X</sub> differed by about 1% (4.70 vs 4.65). Phabs "
          "power-law N<sub>H</sub> was 7.6% higher; Γ = 4.82 was 14.6% steeper (Table 3). At 4 keV, "
          "arfcorr/REEF f<sub>PSF</sub> = 0.1092/0.1077 differed by 1.4%. Rea described Swift/XRT "
          "ARF PSF corrections; ACIS spectra/responses used specextract, but the paper did not state "
          "whether an ACIS aperture correction was applied. ARF scales 0.8/1.0/1.2 gave "
          "L<sub>X</sub> = 5.87/4.70/3.92 × 10³⁵ erg s⁻¹ and R<sub>BB</sub> = 3.35/2.99/2.75 km "
          "(χ²/dof unchanged). Rea’s L<sub>X</sub> = 1.6 × 10³⁵ would require f<sub>PSF</sub> 2.94× "
          "larger (4.70/1.60); neither REEF nor arfcorr supports this, leaving the gap unresolved. "
          "Rea’s Table 1 rate (0.274 ± 0.004 counts s⁻¹) was 3.4× our 0.0808 counts s⁻¹ annular rate. "
          "At 4 keV, REEF f<sub>PSF</sub> = 0.1077 gives 0.0808/0.1077 = 0.7502 ≈ 0.75 counts s⁻¹; "
          "its 0.0778–0.1242 range over 1–10 keV makes this approximate. 0.75/0.274 = 2.74 (about "
          "2.7×), but Rea’s rate region and band were unspecified, so the cause was unknown. Absorbed "
          "phabs 1–10 keV flux was 15.10 × 10⁻¹² erg cm⁻² s⁻¹ (2.4× Rea’s 6.3 × 10⁻¹²); the "
          "retained-band fraction was not measured. Background and N<sub>H</sub>-continuum covariance "
          "remained limitations."),
        p("8  Conclusion", "h1"),
        p("The spectrum was consistent with strongly absorbed thermal emission from a hot neutron-star "
          "region (kT = 0.81 keV, R<sub>BB</sub> = 2.99 km; χ²/dof = 23.06/24), although the similar "
          "power-law fit (20.88/24; Δχ² = 2.18) left the continuum model-dependent. "
          "N<sub>H</sub> ≈ 1.1 × 10²³ cm⁻² indicated heavy Galactic Centre absorption and was "
          "within 7% of Rea et al. L<sub>X</sub> = 4.70 × 10³⁵ erg s⁻¹ (2.9 times Rea et al.) "
          "was 96 times the 4.9 × 10³³ erg s⁻¹ spin-down power, so rotation alone could not power it "
          "and magnetic energy was favored; its 1–10 keV value depended on extrapolation beyond "
          "2.26–6.69 keV. kT was 8% lower, and R<sub>BB</sub> and absorbed flux were 2.1 and 2.4 "
          "times higher. Neither the PSF nor ±20% ARF-scale test resolved the luminosity discrepancy."),
        p("References", "h2"),
        p(f"[1] Rea, N., et al. (2013). <i>ApJL</i>, 775, L34. {href(reference_link, 'doi:10.1088/2041-8205/775/2/L34')}.", "ref"),
        p(f"[2] Chandra X-ray Center. {href('https://cxc.cfa.harvard.edu/cdaftp/byobsid/6/14946/primary/', 'Chandra Data Archive, ObsID 14946')}.", "ref"),
        p(f"[3] Chandra X-ray Center. {href('https://cxc.harvard.edu/csc2/cli/', 'Chandra Source Catalog')}.", "ref"),
        p(f"[4] Chandra X-ray Center. {href('https://cxc.harvard.edu/ciao/threads/pointlike/', 'CIAO point-source spectroscopy guide')}.", "ref"),
        p(f"[5] Chandra X-ray Center. {href('https://cxc.harvard.edu/ciao4.17/ahelp/psf.html', 'CIAO PSF calibration interface')}.", "ref"),
        p(f"[6] NASA HEASARC. {href('https://heasarc.gsfc.nasa.gov/docs/software/xspec/python/html/quick.html', 'XSPEC and PyXspec documentation')}.", "ref"),
        p(f"[7] Chandra X-ray Center. {href('https://cxc.harvard.edu/ciao4.18/ahelp/pileup_map.html', 'CIAO pileup_map help')}.", "ref"),
    ]
    document = SimpleDocTemplate(
        str(pdf_path), pagesize=A4, leftMargin=49, rightMargin=49, topMargin=42, bottomMargin=47,
        title="X-ray spectroscopy of the Galactic Centre magnetar",
        author="SPA660 Project 8 group",
    )
    def first_page(canvas, doc):
        footer(canvas, doc)
        if args.code_url:
            url = escape(args.code_url)
            link = p(f"<b>Project code on GitHub</b>&nbsp;&nbsp; {href(url, url)}", "ref")
            link.wrap(doc.width, 12)
            link.drawOn(canvas, doc.leftMargin + 6, 44)

    document.build(story, onFirstPage=first_page, onLaterPages=footer)


if __name__ == "__main__":
    main()
