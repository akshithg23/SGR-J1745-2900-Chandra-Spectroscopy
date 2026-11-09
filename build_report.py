"""Generate the SPA660 response-fitting report template."""

import json
import math
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


ROOT = Path(__file__).resolve().parent
FONT_PATH = Path("C:/Windows/Fonts")
pdfmetrics.registerFont(TTFont("Report", str(FONT_PATH / "calibri.ttf")))
pdfmetrics.registerFont(TTFont("ReportBold", str(FONT_PATH / "calibrib.ttf")))
pdfmetrics.registerFontFamily("Report", normal="Report", bold="ReportBold",
                              italic="Report", boldItalic="ReportBold")

INK = colors.HexColor("#202833")
MUTED = colors.HexColor("#485668")
ACCENT = colors.HexColor("#235a72")
PALE = colors.HexColor("#e8f1f5")
GRID = colors.HexColor("#d6e0e6")
styles = {
    "title": ParagraphStyle("Title", fontName="ReportBold", fontSize=30, leading=32, textColor=INK, spaceAfter=12),
    "subtitle": ParagraphStyle("Subtitle", fontName="Report", fontSize=14, leading=17, textColor=MUTED, spaceAfter=15),
    "kicker": ParagraphStyle("Kicker", fontName="ReportBold", fontSize=10, leading=12, textColor=ACCENT, spaceAfter=10),
    "h1": ParagraphStyle("Heading", fontName="ReportBold", fontSize=19, leading=23, textColor=INK, spaceBefore=8, spaceAfter=8),
    "h2": ParagraphStyle("Subheading", fontName="ReportBold", fontSize=12, leading=15, textColor=INK, spaceBefore=8, spaceAfter=4),
    "body": ParagraphStyle("Body", fontName="Report", fontSize=11.5, leading=15.5, textColor=INK, spaceAfter=8),
    "caption": ParagraphStyle("Caption", fontName="Report", fontSize=10, leading=13.5, textColor=MUTED, spaceBefore=5, spaceAfter=9),
    "cell": ParagraphStyle("Table", fontName="Report", fontSize=10.5, leading=13, textColor=INK, alignment=TA_LEFT),
    "ref": ParagraphStyle("Reference", fontName="Report", fontSize=10.5, leading=13.5, textColor=INK, spaceAfter=7),
}


def p(value, style="body"):
    return Paragraph(value, styles[style])


def href(url, label):
    return f'<link href="{url}" color="#235a72">{label}</link>'


def table(rows, widths):
    obj = Table([[p(v, "cell") for v in row] for row in rows], colWidths=widths, hAlign="LEFT")
    obj.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), PALE),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f8fa")]),
        ("GRID", (0, 0), (-1, -1), 0.4, GRID),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 9),
        ("RIGHTPADDING", (0, 0), (-1, -1), 9),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    return obj


def graphic(filename, width):
    item = Image(str(ROOT / filename))
    item.drawHeight = item.imageHeight * width / item.imageWidth
    item.drawWidth = width
    return item


def footer(canvas, doc):
    canvas.setFont("Report", 9)
    canvas.setFillColor(MUTED)
    canvas.drawString(49, 31, "SPA660  •  Chandra ObsID 14946")
    canvas.drawRightString(A4[0] - 49, 31, str(doc.page))


def main():
    summary = json.loads((ROOT / "event_summary.json").read_text())
    ns = summary["source_counts_1_10"]
    nb = summary["background_counts_1_10"]
    net = summary["net_counts_1_10"]
    sigma = summary["net_counts_error_1_10"]
    alpha = summary["area_ratio"]
    live = summary["exposure_s"]
    scale = 4 * math.pi * (8.3 * 3.08567758149e21)**2
    benchmark_flux = 1.6e35 / scale
    source_link = "https://cxc.cfa.harvard.edu/cdaftp/byobsid/6/14946/primary/"

    story = [
        p("SPA660 ASTROPHYSICS PROJECT  /  PROJECT 8", "kicker"),
        p("X-ray spectroscopy of the<br/>Galactic Centre magnetar", "title"),
        p("SGR J1745−2900 in Chandra ObsID 14946", "subtitle"),
        p("SPA660 Project 8 group"),
        p("<b>Analysis draft</b> · Bracketed fields require CIAO and XSPEC output", "caption"),
        Spacer(1, 8),
        p("Abstract", "h1"),
        p(f"We examined the public Chandra ACIS-S Level 2 event list for SGR J1745−2900. A 1–8 keV "
          f"image resolves the magnetar beside Sgr A*. After excluding the central 0.75 arcsec and "
          f"scaling the local background by the extraction area, the 1–10 keV annulus yielded "
          f"{net:,.0f} ± {sigma:.0f} net counts. The detector-energy distribution peaks in the "
          "3.50–3.75 keV bin. The response-folded comparison gives [preferred model or comparable fits = ?], "
          "with [χ²/dof = ?], [N<sub>H</sub> = ? cm⁻²] and [kT = ? keV or Γ = ?]. "
          "The intrinsic 1–10 keV luminosity at 8.3 kpc is [L<sub>X</sub> = ? erg s⁻¹]."),
        p("1  Target and observation", "h1"),
        p("SGR J1745−2900 is a transient magnetar near the Galactic Centre. Its published "
          "position lies 2.4 ± 0.3 arcsec from Sgr A* [1]. The proximity of the two sources "
          "and the diffuse emission in the field make a small extraction aperture and local "
          "background estimate important."),
        table([
            ["<b>Observation property</b>", "<b>Value</b>"],
            ["Chandra observation", "ObsID 14946 · 2 July 2013"],
            ["Detector and mode", "ACIS-S3 (CCD 7) · FAINT · 1/8 subarray"],
            ["Exposure and frame time", f"{live/1000:.3f} ks · {summary['frame_time_s']:.5f} s"],
            ["Archive calibration", f"CALDB {summary['caldb']} · Level 2 event list [2]"],
        ], [192, 305]),
        Spacer(1, 12),
        p("2  Event selection", "h1"),
        p("We selected CCD 7 events and used the celestial coordinate transform in the FITS "
          "event table to register the image. The source aperture is a 0.75–2.0 arcsec annulus "
          "centred at RA 266.4173708° and Dec −29.0082889° (J2000), following the geometry "
          "of Rea et al. [1]. Excluding the bright core reduces pileup, in which multiple "
          "photons within a frame can be recorded as one event or rejected."),
        p("We estimated the local background in a concentric 2.5–8.0 arcsec annulus, excluding "
          "Sgr A* and nearby catalogued sources [3]. The background region samples the local "
          "diffuse emission. Wings from nearby sources can still contribute to the source aperture."),
        PageBreak(),
        p("3  Galactic Centre image", "h1"),
        graphic("field.png", 492),
        p("<b>Figure 1.</b> 1–8 keV archive event image. Solid cyan circles bound the source annulus. "
          "Dashed yellow circles bound the background annulus and dotted white circles mark its "
          "point-source exclusions. The axes are relative to the published magnetar position. "
          "The colour stretch saturates bright central pixels for visibility.", "caption"),
        p("The bright core falls at the magnetar’s published position. The projected separation "
          "of the adopted magnetar and Sgr A* coordinates is "
          f"{summary['separation_from_reference_sgra_arcsec']:.2f} arcsec. This separation is "
          "calculated from reference coordinates [1,3]."),
        p("Aperture photometry", "h2"),
        table([
            ["<b>Quantity</b>", "<b>1–10 keV result</b>"],
            ["Source annulus events", f"{ns:,}"],
            ["Masked background annulus events", f"{nb:,}"],
            ["Source-to-background area ratio", f"{alpha:.5f}"],
            ["Net source counts", f"{net:,.0f} ± {sigma:.0f}"],
            ["Net annular count rate", f"{net/live:.4f} ± {sigma/live:.4f} counts s⁻¹"],
        ], [286, 211]),
        Spacer(1, 8),
        p("The uncertainty is the propagated Poisson counting error. We computed net counts as "
          f"N<sub>s</sub> − αN<sub>b</sub>, with α = {alpha:.5f} from the unmasked geometric areas. "
          "The rate describes the aperture after background subtraction and before correction for "
          "photons outside the annulus.", "caption"),
        PageBreak(),
        p("4  Detector spectrum and aperture correction", "h1"),
        graphic("counts.png", 497),
        p("<b>Figure 2.</b> Background-subtracted count-rate density from 1–10 keV. "
          "Counts and their propagated errors were divided by the 18.202 ks exposure and "
          "the 0.25 keV bin width. The orange histogram shows the area-scaled background. "
          "The x-axis uses the archive event ENERGY column.", "caption"),
        p("The signal is strongest between about 3 and 5 keV and declines sharply above "
          "6 keV. The observed shape depends on the source continuum, line-of-sight absorption, "
          "the ACIS collecting area and energy redistribution. An instrument response is required "
          "to infer the intrinsic spectrum."),
        p("CIAO reprocessing and pileup", "h2"),
        p("Use the complete ObsID 14946 archive dataset as input to <b>chandra_repro</b> "
          "[CIAO version = ?; reprocessing CALDB = ?]. Confirm the adopted position in DS9 using the "
          "1–8 keV CCD 7 image. A full-band native-pixel CCD 7 image is the input to "
          "<b>pileup_map</b>. The maxima are [core = ? counts frame⁻¹] and [annulus = ? counts frame⁻¹]. "
          "Rea et al. used a 0.75 arcsec inner cut to reduce core pileup [1]. The 2.0 arcsec outer "
          "radius limits contamination from Sgr A* and diffuse emission."),
        p("Spectrum and point spread function", "h2"),
        p("CIAO <b>specextract</b> is configured to produce source and background PI spectra, "
          "the RMF and the ARF. The energy-dependent <b>arfcorr</b> factor multiplies the effective "
          "area by the PSF fraction in the annulus. At 4 keV, [f<sub>PSF</sub> = ?] and "
          "[1/f<sub>PSF</sub> = ?]. Across 1–10 keV the fraction spans [range = ?]. "
          "Use the corrected ARF in XSPEC to account for the missing core and wings."),
        PageBreak(),
        p("5  XSPEC fits and reference comparison", "h1"),
        p("Fit the source groups with at least 50 counts over 1–10 keV using standard χ² weighting. "
          "Use <b>tbabs*bbody</b> and <b>tbabs*powerlaw</b> with a free column and with "
          "N<sub>H</sub> fixed at 10.3 × 10²² cm⁻². XSPEC uses <b>angr</b> abundances. "
          "<b>tbabs</b> uses its own cross-sections. The <b>phabs</b> cross-check uses "
          "<b>bcmc</b> for comparison with the paper [1,5]. Parameter intervals use Δχ² = 1 "
          "for one parameter. Flux intervals come from <b>cflux</b> around the full model "
          "and around the continuum alone for absorbed and intrinsic flux, respectively."),
        table([
            ["<b>Continuum and column</b>", "<b>N<sub>H</sub> (10²² cm⁻²)</b>", "<b>kT (keV) or Γ</b>", "<b>χ²/dof</b>"],
            ["Blackbody, free", "[N<sub>H</sub> = ?]", "[kT = ?]", "[? / ?]"],
            ["Blackbody, frozen", "10.3 (fixed)", "[kT = ?]", "[? / ?]"],
            ["Power law, free", "[N<sub>H</sub> = ?]", "[Γ = ?]", "[? / ?]"],
            ["Power law, frozen", "10.3 (fixed)", "[Γ = ?]", "[? / ?]"],
        ], [156, 114, 117, 110]),
        Spacer(1, 7),
        p("Table 2. Report best-fit values with 68% intervals for the free parameters. Compare "
          "the two continua at the same column setting using χ² and residual structure.", "caption"),
        p("Comparison with Rea et al.", "h2"),
        table([
            ["<b>Quantity and unit</b>", "<b>Free-column fit</b>", "<b>Rea et al.</b>", "<b>Δ (%)</b>"],
            ["Blackbody N<sub>H</sub> (10²² cm⁻²)", "[?]", "10.3 ± 0.4", "[?]"],
            ["Power-law N<sub>H</sub> (10²² cm⁻²)", "[?]", "17.3 ± 0.4", "[?]"],
            ["Blackbody kT (keV)", "[?]", "0.88 ± 0.02", "[?]"],
            ["Power-law Γ", "[?]", "4.2 ± 0.1", "[?]"],
            ["Absorbed flux (10⁻¹² erg cm⁻² s⁻¹)", "[?]", "6.3 ± 0.2", "[?]"],
            ["Intrinsic flux (10⁻¹¹ erg cm⁻² s⁻¹)", "[?]", f"≈ {benchmark_flux/1e-11:.1f} (derived)", "[?]"],
            ["Luminosity (10³⁵ erg s⁻¹)", "[?]", "1.6", "[?]"],
        ], [214, 100, 109, 74]),
        p("Table 3. Flux and luminosity rows use the blackbody benchmark for ObsID 14946. "
          "The reference intrinsic flux is derived from L/(4πd²). The published columns summarize "
          "three observations. Δ = 100 × (fit − reference)/reference.", "caption"),
        PageBreak(),
        p("6  Spectral results and physical interpretation", "h1"),
        p("<b>Continuum.</b> With N<sub>H</sub> free, the difference "
          "χ²<sub>BB</sub> − χ²<sub>PL</sub> is [Δχ² = ?] at equal degrees of freedom. "
          "The residuals show [structure and relative fit quality = ?]. Freezing the column gives "
          "[change in fit quality and continuum parameters = ?]. The blackbody fit gives "
          "[kT = ? keV] and [apparent R<sub>BB</sub> = ? km]. If the fit is adequate and "
          "the radius is small compared with a neutron star, it supports emission from a localized hot region."),
        p("<b>Absorption.</b> The fitted [model = ?] column of [N<sub>H</sub> = ? cm⁻²] "
          "implies [strength of soft X-ray attenuation = ?] along the Galactic Centre sight line. "
          "The model gives [absorbed flux = ?] and [intrinsic flux = ?] in erg cm⁻² s⁻¹ over 1–10 keV."),
        p("<b>Luminosity and energy distribution.</b> The intrinsic luminosity at 8.3 kpc is "
          "[L<sub>X</sub> = ? erg s⁻¹], or [ratio = ?] times the published spin-down power of "
          "4.9 × 10³³ erg s⁻¹ [1]. If the luminosity exceeds this power, rotation alone "
          "cannot sustain the emission. The absorbed and intrinsic energy flux densities at 4 keV are "
          "[F<sub>E</sub> = ? and ? erg cm⁻² s⁻¹ keV⁻¹]. The model F<sub>E</sub> maximum "
          "in 1–10 keV is [energy = ? keV; interior or band edge = ?]."),
        p("<b>Agreement with the reference.</b> The measured [kT, Γ, column or flux = ?] "
          "differs by [Δ = ?%] and is [consistent or discrepant given the intervals = ?]. "
          "The fit using <b>phabs</b> gives [change relative to <b>tbabs</b> = ?]."),
        p("<b>Figure 3.</b> [Insert selected <b>*_fit.png</b>.] The upper panel shows the "
          "grouped count spectrum and response-folded model. The lower panel shows residuals in σ.", "caption"),
        p("<b>Figure 4.</b> [Insert corresponding <b>*_flux_density.png</b>.] "
          "The unfolded points depend on the fit. The model curves show absorbed and intrinsic F<sub>E</sub>.", "caption"),
        p("Limitations", "h2"),
        p("A single epoch constrains the continuum and absorption less tightly than the "
          "three-observation fit. The absorption column covaries with kT or Γ. The narrow PSF "
          "aperture, diffuse background, nearby-source wings and residual pileup can shift the "
          "normalization or shape. The two continua are non-nested and a small χ² difference "
          "alone gives limited evidence for a preference."),
        p("7  Conclusion", "h1"),
        p(f"The archive observation resolves SGR J1745−2900 and yields {net:,.0f} ± {sigma:.0f} "
          "net counts in the pileup-excluding annulus. The aperture-corrected spectrum is "
          "[fit conclusion = ?], with [N<sub>H</sub> = ? cm⁻²], [kT = ? keV or Γ = ?] "
          "and [L<sub>X</sub> = ? erg s⁻¹]. The main agreement or difference with Rea et al. is "
          "[comparison = ?]."),
        p("References", "h1"),
        p("[1] Rea, N., et al. (2013). A strongly magnetized pulsar within grasp of the Milky "
          f"Way’s supermassive black hole. <i>ApJL</i>, 775, L34. "
          f"{href('https://doi.org/10.1088/2041-8205/775/2/L34', 'doi:10.1088/2041-8205/775/2/L34')}.", "ref"),
        p(f"[2] Chandra X-ray Center. {href(source_link, 'Chandra Data Archive, ObsID 14946')}.", "ref"),
        p(f"[3] Chandra X-ray Center. {href('https://cxc.harvard.edu/csc2/cli/', 'Chandra Source Catalog')}. "
          "Source positions within 10 arcsec of SGR J1745−2900.", "ref"),
        p(f"[4] Chandra X-ray Center. {href('https://cxc.harvard.edu/ciao/threads/pointlike/', 'CIAO point-source spectroscopy guide')}.", "ref"),
        p(f"[5] NASA HEASARC. {href('https://heasarc.gsfc.nasa.gov/docs/software/xspec/python/html/quick.html', 'XSPEC and PyXspec documentation')}.", "ref"),
    ]
    output = Path(sys.argv[1])
    output.parent.mkdir(parents=True, exist_ok=True)
    document = SimpleDocTemplate(str(output), pagesize=A4, leftMargin=49, rightMargin=49,
                                 topMargin=42, bottomMargin=47,
                                 title="X-ray spectroscopy of the Galactic Centre magnetar",
                                 author="SPA660 Project 8 group")
    document.build(story, onFirstPage=footer, onLaterPages=footer)


if __name__ == "__main__":
    main()
