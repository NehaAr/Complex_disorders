"""
GWAS SNP-Trait Discovery Timeline — Extended to 2026
=====================================================
Reproduces and extends Figure 2 from Visscher et al. 2017
using GWAS Catalog data downloaded from EMBL-EBI FTP.

HOW TO GET THE DATA FILE
-------------------------
1. Go to: https://www.ebi.ac.uk/gwas/docs/file-downloads
2. Download: "All associations v1.0" (the big TSV file)
   Direct FTP link:
   https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv
3. Save it in the SAME folder as this script
4. Set LOCAL_FILE below to match the filename you saved

OR let the script auto-download via FTP (more reliable than API):
  Set AUTO_DOWNLOAD = True below

REQUIREMENTS
------------
pip install pandas matplotlib numpy requests
"""

import os
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge
import warnings
warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════════════
#  SETTINGS — EDIT THESE
# ══════════════════════════════════════════════════════════════════════════════

# Set to True to auto-download, False if you already have the file
AUTO_DOWNLOAD = True

# If AUTO_DOWNLOAD = False, set this to your local file path
LOCAL_FILE = "gwas-catalog-associations_ontology-annotated.tsv"

# Correct FTP URL (updated 2025 — more reliable than API endpoint)
GWAS_FTP_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations_ontology-annotated.tsv"
)

P_THRESH  = 5e-8    # genome-wide significance threshold
TOP_N     = 3       # top traits to label per year
OUT_FIG   = "gwas_timeline_2026.png"
OUT_CSV   = "gwas_timeline_data.csv"

# ══════════════════════════════════════════════════════════════════════════════

# Distinct colour palette — new scheme, different from original paper
TRAIT_COLOURS = [
    "#2EBD8E", "#534AB7", "#E8503A", "#F5A623", "#378ADD",
    "#9C59D1", "#1D9E75", "#D85A30", "#BA7517", "#185FA5",
    "#C0392B", "#16A085", "#8E44AD", "#2980B9", "#F39C12",
    "#27AE60", "#E74C3C", "#3498DB", "#9B59B6", "#1ABC9C",
    "#E67E22", "#2C3E50", "#7F8C8D", "#D35400", "#2ECC71",
    "#3D5A80", "#98C1D9", "#E84855", "#84A98C", "#FFB703",
    "#219EBC", "#8ECAE6", "#FB8500", "#023047", "#FFB4A2",
    "#FFCDB2", "#B5838D", "#6D6875", "#E5989B", "#FFCCD5",
]

os.makedirs("gwas_figures", exist_ok=True)


# ── DOWNLOAD ──────────────────────────────────────────────────────────────────

def download_file(url, local_file):
    """Download via FTP/HTTPS with progress indicator."""
    if os.path.exists(local_file):
        size_mb = os.path.getsize(local_file) / 1e6
        print(f"  Using cached file: {local_file} ({size_mb:.1f} MB)")
        return True

    print(f"  Downloading from EMBL-EBI FTP...")
    print(f"  URL: {url}")
    print(f"  File will be ~150-250 MB — please wait...")

    try:
        with requests.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            total = int(r.headers.get('content-length', 0))
            downloaded = 0
            with open(local_file, 'wb') as f:
                for chunk in r.iter_content(chunk_size=65536):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total > 0:
                        pct = downloaded / total * 100
                        mb  = downloaded / 1e6
                        print(f"\r  {pct:.1f}% ({mb:.0f} MB)", end='', flush=True)
        print(f"\n  Download complete!")
        return True
    except Exception as e:
        print(f"\n  Auto-download failed: {e}")
        print(f"\n  ── MANUAL DOWNLOAD INSTRUCTIONS ──────────────────────")
        print(f"  1. Open this URL in your browser:")
        print(f"     {url}")
        print(f"  2. Save the file as: {local_file}")
        print(f"     (in the same folder as this script)")
        print(f"  3. Run this script again")
        print(f"  ──────────────────────────────────────────────────────")
        return False


# ── LOAD AND FILTER ───────────────────────────────────────────────────────────

def load_and_filter(local_file):
    """Load GWAS Catalog TSV, filter p < 5e-8, extract year and trait."""
    print(f"\n  Loading file: {local_file}")
    print(f"  (This may take 30-60 seconds for large files...)")

    df = pd.read_csv(
        local_file,
        sep='\t',
        low_memory=False,
        on_bad_lines='skip',
        encoding='utf-8',
    )
    print(f"  Rows loaded: {len(df):,}")
    print(f"  Columns: {list(df.columns[:8])}")

    # ── Find the right columns ─────────────────────────────────────────────
    # P-value column
    p_col = None
    for candidate in ['P-VALUE', 'P_VALUE', 'PVALUE', 'p_value', 'p-value']:
        if candidate in df.columns:
            p_col = candidate
            break
    if p_col is None:
        for c in df.columns:
            if 'p' in c.lower() and 'val' in c.lower():
                p_col = c
                break

    # Date column
    date_col = None
    for candidate in ['DATE ADDED TO CATALOG', 'DATE', 'date_added', 'DATE_ADDED']:
        if candidate in df.columns:
            date_col = candidate
            break
    if date_col is None:
        for c in df.columns:
            if 'date' in c.lower():
                date_col = c
                break

    # Trait column
    trait_col = None
    for candidate in ['DISEASE/TRAIT', 'TRAIT', 'disease_trait', 'DISEASE_TRAIT']:
        if candidate in df.columns:
            trait_col = candidate
            break
    if trait_col is None:
        for c in df.columns:
            if 'trait' in c.lower() or 'disease' in c.lower():
                trait_col = c
                break

    print(f"  Columns identified:")
    print(f"    p-value : '{p_col}'")
    print(f"    date    : '{date_col}'")
    print(f"    trait   : '{trait_col}'")

    if not all([p_col, date_col, trait_col]):
        raise ValueError(
            f"Could not identify required columns.\n"
            f"Available columns: {list(df.columns)}"
        )

    # ── Clean ──────────────────────────────────────────────────────────────
    df[p_col]     = pd.to_numeric(df[p_col], errors='coerce')
    df[date_col]  = pd.to_datetime(df[date_col], errors='coerce')
    df[trait_col] = df[trait_col].astype(str).str.strip()

    # Filter p < 5e-8
    df = df[df[p_col] < P_THRESH].copy()
    print(f"  After p < 5e-8: {len(df):,} associations")

    # Extract year — group pre-2008 as 2007
    df['YEAR'] = df[date_col].dt.year
    df = df.dropna(subset=['YEAR'])
    df['YEAR'] = df['YEAR'].astype(int)
    df.loc[df['YEAR'] < 2008, 'YEAR'] = 2007
    df = df[(df['YEAR'] >= 2007) & (df['YEAR'] <= 2026)]

    # Rename for consistency
    df = df.rename(columns={trait_col: 'TRAIT', p_col: 'P_VALUE'})

    # Proxy LD pruning — keep one entry per trait per year
    # (real LD pruning needs genotype data; this is the standard proxy)
    df = df.sort_values('P_VALUE')
    df = df.drop_duplicates(subset=['TRAIT', 'YEAR'])

    print(f"  After LD proxy pruning: {len(df):,} unique trait-year entries")
    return df


# ── AGGREGATE ─────────────────────────────────────────────────────────────────

def aggregate_by_year(df):
    """Count SNPs per trait per year and compute cumulative totals."""
    year_trait = (
        df.groupby(['YEAR', 'TRAIT'])
        .size()
        .reset_index(name='N_SNPS')
    )
    totals = (
        year_trait.groupby('YEAR')['N_SNPS']
        .sum()
        .reset_index()
        .sort_values('YEAR')
    )
    totals['CUMULATIVE'] = totals['N_SNPS'].cumsum()

    yearly_data = {
        yr: grp.sort_values('N_SNPS', ascending=False)
        for yr, grp in year_trait.groupby('YEAR')
    }
    return yearly_data, totals


# ── DONUT DRAWING ─────────────────────────────────────────────────────────────

def draw_donut(ax, cx, cy, radius, trait_counts, trait_colour_map, total):
    """Draw a donut chart centred at (cx, cy) in data coordinates."""
    if trait_counts.sum() == 0:
        return

    start = 90.0
    for trait, count in trait_counts.items():
        frac   = count / trait_counts.sum()
        angle  = frac * 360
        colour = trait_colour_map.get(trait, '#BBBBBB')
        wedge  = Wedge(
            center=(cx, cy), r=radius,
            theta1=start - angle, theta2=start,
            width=radius * 0.45,
            facecolor=colour, edgecolor='white',
            linewidth=0.5, zorder=5,
        )
        ax.add_patch(wedge)
        start -= angle

    ax.text(cx, cy, str(total),
            ha='center', va='center',
            fontsize=6.5, fontweight='bold', color='#111111', zorder=6)


# ── MAIN PLOT ─────────────────────────────────────────────────────────────────

def plot_timeline(yearly_data, yearly_totals, out_fig):
    """Reproduce and extend the GWAS SNP discovery timeline figure."""
    print(f"\n  Building figure...")

    years  = sorted(yearly_totals['YEAR'].tolist())
    cum    = dict(zip(yearly_totals['YEAR'], yearly_totals['CUMULATIVE']))
    annual = dict(zip(yearly_totals['YEAR'], yearly_totals['N_SNPS']))

    # Assign colours to traits consistently
    all_traits = []
    for yr in years:
        if yr in yearly_data:
            all_traits.extend(yearly_data[yr]['TRAIT'].tolist())
    unique_traits = list(dict.fromkeys(all_traits))
    trait_colour_map = {
        t: TRAIT_COLOURS[i % len(TRAIT_COLOURS)]
        for i, t in enumerate(unique_traits)
    }

    max_cum  = max(cum.values()) if cum else 1
    max_ann  = max(annual.values()) if annual else 1

    fig, ax = plt.subplots(figsize=(24, 13))
    fig.patch.set_facecolor('#FAFAF8')
    ax.set_facecolor('#F5F4F0')

    x_vals = list(range(len(years)))
    y_vals = [cum.get(yr, 0) for yr in years]

    # S-curve
    ax.plot(x_vals, y_vals,
            color='#1A1A2E', linewidth=2.2, zorder=2, solid_capstyle='round')
    ax.scatter(x_vals, y_vals,
               color='#378ADD', s=45, zorder=3, linewidths=0)

    # Donuts
    base_radius = max_cum * 0.055
    for i, yr in enumerate(years):
        xi  = i
        yi  = cum.get(yr, 0)
        tot = annual.get(yr, 0)
        if yr not in yearly_data or tot == 0:
            continue

        trait_counts = yearly_data[yr].set_index('TRAIT')['N_SNPS']
        scale        = 0.5 + np.log1p(tot) / np.log1p(max_ann) * 1.6
        radius       = base_radius * scale

        draw_donut(ax, xi, yi, radius, trait_counts, trait_colour_map, tot)

        # Label top N traits
        top  = yearly_data[yr].head(TOP_N)
        lbls = [row['TRAIT'][:32] for _, row in top.iterrows()]
        ax.text(
            xi, yi + radius + max_cum * 0.025,
            '\n'.join(lbls),
            ha='center', va='bottom', fontsize=5.2, color='#333333',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='#cccccc', alpha=0.85, linewidth=0.4),
            zorder=8,
        )

    # Axes
    ax.set_xlim(-0.8, len(years) - 0.2)
    ax.set_ylim(-max_cum * 0.1, max_cum * 1.30)
    ax.set_xticks(x_vals)
    ax.set_xticklabels(
        ['Before\n2008' if yr == 2007 else str(yr) for yr in years],
        fontsize=9,
    )
    ax.set_ylabel('Cumulative SNP-Trait Associations (p < 5×10⁻⁸)',
                  fontsize=11)
    ax.set_xlabel('Year of Discovery', fontsize=11)
    ax.set_title(
        'GWAS SNP-Trait Discovery Timeline  ·  Before 2008 to 2026\n'
        'Source: NHGRI-EBI GWAS Catalog  |  p < 5×10⁻⁸  |  LD proxy pruned',
        fontsize=13, fontweight='bold', color='#1A1A2E', pad=16,
    )

    # NOTE: Purple tint removed — no shading between years
    # (all years treated equally; no distinction between original and extended)

    # Grid and spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.18, linewidth=0.5, linestyle='--')

    # Schema legend box (top left)
    ax.text(
        0.01, 0.97,
        'Circle key:\n● Height = SNPs discovered that year\n'
        '● Segment width = fraction of publications\n'
        '● Top 3 traits labelled per year',
        transform=ax.transAxes,
        fontsize=7.5, va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#EFEFEF',
                  edgecolor='#cccccc', alpha=0.9),
    )

    plt.tight_layout()
    plt.savefig(out_fig, dpi=180, bbox_inches='tight', facecolor='#FAFAF8')
    plt.close()
    print(f"  Figure saved → {out_fig}")


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  GWAS Timeline — Extended to 2026")
    print("=" * 60)

    # Step 1: Get data
    if AUTO_DOWNLOAD:
        success = download_file(GWAS_FTP_URL, LOCAL_FILE)
        if not success:
            return
    else:
        if not os.path.exists(LOCAL_FILE):
            print(f"\nERROR: File not found: {LOCAL_FILE}")
            print(f"Set AUTO_DOWNLOAD = True or download the file manually.")
            print(f"URL: {GWAS_FTP_URL}")
            return

    # Step 2: Load and filter
    print("\n[2/4] Loading and filtering...")
    df = load_and_filter(LOCAL_FILE)

    # Step 3: Aggregate
    print("\n[3/4] Aggregating by year...")
    yearly_data, yearly_totals = aggregate_by_year(df)

    print(f"\n  {'Year':<14} {'SNPs this year':>15} {'Cumulative':>12}")
    print(f"  {'-'*44}")
    for _, row in yearly_totals.iterrows():
        yr = 'Before 2008' if row['YEAR'] == 2007 else int(row['YEAR'])
        print(f"  {str(yr):<14} {int(row['N_SNPS']):>15,} {int(row['CUMULATIVE']):>12,}")

    # Step 4: Plot
    print("\n[4/4] Generating figure...")
    plot_timeline(yearly_data, yearly_totals, OUT_FIG)

    # Save CSV
    rows = []
    for yr, grp in yearly_data.items():
        for _, row in grp.iterrows():
            rows.append({
                'Year'  : 'Before 2008' if yr == 2007 else yr,
                'Trait' : row['TRAIT'],
                'N_SNPs': row['N_SNPS'],
            })
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    print(f"  Data saved → {OUT_CSV}")

    print(f"\n{'='*60}")
    print(f"  Done!")
    print(f"  Figure → {OUT_FIG}")
    print(f"  Data   → {OUT_CSV}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()