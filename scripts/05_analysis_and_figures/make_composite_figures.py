#!/usr/bin/env python3
"""Assembles the lettered multi-panel manuscript figures (Fig2, Fig3, Fig6) from the figure scripts' PNGs.

Usage:
    python make_composite_figures.py --bootstrap-figs bootstrap_figs_out --intra-inter intra_vs_inter_out --treatment-phase treatment_phase_out --out figures
"""
import argparse
import os
from PIL import Image, ImageDraw, ImageFont

TARGET_H = 1900   # common panel height before label padding
LABEL_PAD = 90    # space above panels for the letter
GAP = 30


def load_scaled(path):
    im = Image.open(path).convert("RGB")
    w, h = im.size
    return im.resize((int(w * TARGET_H / h), TARGET_H), Image.LANCZOS)


def try_font(size):
    for fp in ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
               "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
        if os.path.exists(fp):
            return ImageFont.truetype(fp, size)
    return ImageFont.load_default()


def build_row(panel_paths, labels, out_path):
    missing = [p for p in panel_paths if not os.path.exists(p)]
    if missing:
        print("[skip] %s -- missing inputs: %s" % (os.path.basename(out_path), missing))
        return
    imgs = [load_scaled(p) for p in panel_paths]
    total_w = sum(im.width for im in imgs) + GAP * (len(imgs) - 1)
    canvas = Image.new("RGB", (total_w, TARGET_H + LABEL_PAD), "white")
    draw = ImageDraw.Draw(canvas)
    font = try_font(70)
    x = 0
    for im, label in zip(imgs, labels):
        canvas.paste(im, (x, LABEL_PAD))
        draw.text((x + 10, 5), label, fill="black", font=font)
        x += im.width + GAP
    canvas.save(out_path, dpi=(300, 300))
    print("Wrote", out_path, canvas.size)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bootstrap-figs", default="bootstrap_figs_out")
    ap.add_argument("--intra-inter", default="intra_vs_inter_out")
    ap.add_argument("--treatment-phase", default="treatment_phase_out")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    j = os.path.join
    build_row([j(a.bootstrap_figs, "fig1_bootstrap_support/support_all_patients.png"),
               j(a.bootstrap_figs, "fig2_clone_normal_jaccard/jaccard_all_patients.png"),
               j(a.bootstrap_figs, "fig3_intra_inter_hamming/hamming_all_patients.png")],
              ["a", "b", "c"], j(a.out, "Fig2.png"))
    build_row([j(a.intra_inter, "violin_all_patients.png"),
               j(a.intra_inter, "stddev_boxplot_by_cancer_type.png"),
               j(a.bootstrap_figs, "fig4_clone_size_divergence/divergence_all_patients.png")],
              ["a", "b", "c"], j(a.out, "Fig3.png"))
    build_row([j(a.treatment_phase, "wasserstein_within_vs_cross.png"),
               j(a.treatment_phase, "gene_degree_shift/top_gene_degree_shift.png")],
              ["a", "b"], j(a.out, "Fig6.png"))


if __name__ == "__main__":
    main()
