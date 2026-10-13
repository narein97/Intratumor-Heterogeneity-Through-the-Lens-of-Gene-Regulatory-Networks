"""Builds the per-cancer-type supplementary composites (Supplementary Figs S11-S18).

Usage:
    python make_supplementary_composites.py --bootstrap-figs bootstrap_figs_out --out supp_figs
"""
import os
from PIL import Image, ImageDraw, ImageFont
"""Build the per-cancer-type supplementary composites (Supplementary Figs S11-S18) from the per-cancer-type PNGs written by
bootstrap_figures_efficient.py and intra_vs_inter_efficient.py.

Usage: python make_supplementary_composites.py --bootstrap-figs bootstrap_figs_out --out supp_figs
Writes SuppExt_{prop_1..3,count_1..3,bootstrap_support,hamming}.png (panel letters added; 4 cancer types per clone-size figure, 12 panels
for the bootstrap-support and Hamming figures). 'Ovar' is the directory name used for ovarian in the cohort-wide outputs.
"""
import argparse
ap=argparse.ArgumentParser()
ap.add_argument("--bootstrap-figs",required=True,help="bootstrap_figures_efficient.py output directory")
ap.add_argument("--out",default="supp_figs")
_a=ap.parse_args()
B=_a.bootstrap_figs; OUT=_a.out
os.makedirs(OUT,exist_ok=True)
CT=["Brain","Breast","Colorectal","Head_and_Neck","Hematologic","Kidney","Lung","Neuroendocrine","Ovar","Pancreas","Prostate","Skin"]
def font(s):
    for fp in ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf","/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
        if os.path.exists(fp): return ImageFont.truetype(fp,s)
    return ImageFont.load_default()
LET="abcdefghijkl"
def stack(paths,letters,out,fs=170,pad=230,gap=40):
    ims=[Image.open(p).convert("RGB") for p in paths]
    W=max(i.width for i in ims); H=sum(i.height+pad for i in ims)+gap*(len(ims)-1)
    c=Image.new("RGB",(W,H),"white"); d=ImageDraw.Draw(c); f=font(fs); y=0
    for im,l in zip(ims,letters):
        d.text((20,y+10),l,fill="black",font=f); c.paste(im,(0,y+pad)); y+=im.height+pad+gap
    c.save(out,dpi=(300,300)); print(out,c.size)
def grid(paths,letters,out,cols,th,fs=130,pad=170,gap=30):
    ims=[Image.open(p).convert("RGB") for p in paths]
    ims=[i.resize((int(i.width*th/i.height),th),Image.LANCZOS) for i in ims]
    cw=max(i.width for i in ims); rows=(len(ims)+cols-1)//cols
    c=Image.new("RGB",(cols*cw+gap*(cols-1),rows*(th+pad)+gap*(rows-1)),"white"); d=ImageDraw.Draw(c); f=font(fs)
    for k,(im,l) in enumerate(zip(ims,letters)):
        r,q=divmod(k,cols); x=q*(cw+gap); y=r*(th+pad+gap)
        d.text((x+10,y+5),l,fill="black",font=f); c.paste(im,(x+(cw-im.width)//2,y+pad))
    c.save(out,dpi=(300,300)); print(out,c.size)
for kind,pat in [("prop","fig4_clone_size_divergence/clone_size_proportion_%s.png"),("count","fig4_clone_size_divergence/clone_size_count_%s.png")]:
    for k in range(3):
        cs=CT[4*k:4*k+4]
        stack([os.path.join(B,pat%c) for c in cs],list(LET[4*k:4*k+4]),os.path.join(OUT,f"SuppExt_{kind}_{k+1}.png"))
grid([os.path.join(B,"fig1_bootstrap_support/support_%s.png"%c) for c in CT],list(LET),os.path.join(OUT,"SuppExt_bootstrap_support.png"),3,1500)
grid([os.path.join(B,"fig3_intra_inter_hamming/hamming_%s.png"%c) for c in CT],list(LET),os.path.join(OUT,"SuppExt_hamming.png"),3,1500)
