#!/usr/bin/env python3
import re, math, argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PAGE_W = 26  # inches (wide landscape)
PAGE_H = 12

def _adjust_layout(fig, ax=None):
    fig.subplots_adjust(left=0.06, right=0.96, top=0.92, bottom=0.09)
    if ax is not None:
        ax.set_position([0.055, 0.12, 0.89, 0.78])

def _page_number(fig, num):
    fig.text(0.5, 0.035, f"{num}", ha="center", va="center", fontsize=9)

def wrap_text(s, width):
    if s is None:
        return ""
    return "\n".join([str(s)[i:i+width] for i in range(0, len(str(s)), width)]) if width else str(s)

def draw_table(ax, df: pd.DataFrame, title: str, col_wrap: dict = None, fontsize=8, scale=(1.5, 1.5), header_wrap=14):
    ax.axis('off'); ax.set_title(title, fontsize=12, pad=10)
    if col_wrap:
        df = df.copy()
        for c,w in col_wrap.items():
            if c in df.columns:
                df[c] = df[c].map(lambda x: wrap_text(x, w))
    headers = [wrap_text(h, header_wrap) for h in df.columns]
    tbl = ax.table(cellText=df.values, colLabels=headers, loc='center', cellLoc='center', colLoc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(fontsize); tbl.scale(*scale)
    for (r,c), cell in tbl.get_celld().items():
        cell.set_edgecolor('0.6'); cell.set_linewidth(0.4)
        if r == 0:
            fp = cell.get_text().get_fontproperties(); fp.set_weight('bold'); cell.get_text().set_fontproperties(fp)

def render_wide_table(pdf, title, df, cols_per_page=15, col_wrap=None, fontsize=7, scale=(1.45,1.45), page_no_ref=[1]):
    if df.empty:
        fig, ax = plt.subplots(figsize=(PAGE_W, PAGE_H)); _adjust_layout(fig, ax); ax.axis('off')
        ax.text(0.5, 0.5, f"{title}: no data", ha='center', va='center')
        _page_number(fig, page_no_ref[0]); pdf.savefig(fig); plt.close(fig); page_no_ref[0] += 1
        return
    cols = list(df.columns)
    for i in range(0, len(cols), cols_per_page):
        sub = cols[i:i+cols_per_page]
        fig, ax = plt.subplots(figsize=(PAGE_W, PAGE_H)); _adjust_layout(fig, ax)
        draw_table(ax, df[sub], f"{title} (cols {i+1}-{i+len(sub)}/{len(cols)})", col_wrap=col_wrap, fontsize=fontsize, scale=scale, header_wrap=14)
        _page_number(fig, page_no_ref[0]); pdf.savefig(fig); plt.close(fig); page_no_ref[0] += 1

def collect_files(base: Path, prefix: str):
    def safe_sort(glob_pat, regex_pat):
        files = list(base.glob(glob_pat))
        def keyfun(p):
            m = re.search(regex_pat, p.name)
            return int(m.group(1)) if m else 10**9
        return sorted(files, key=keyfun)
    srv = safe_sort(f"{prefix}_*_codon_srv.out.txt", rf"{re.escape(prefix)}_(\d+)_codon_srv\.out\.txt")
    mult2 = safe_sort(f"{prefix}_*_codon_srv_multDoub.out.txt", rf"{re.escape(prefix)}_(\d+)_codon_srv_multDoub\.out\.txt")
    mult23 = safe_sort(f"{prefix}_*_codon_srv_multDoubTrip.out.txt", rf"{re.escape(prefix)}_(\d+)_codon_srv_multDoubTrip\.out\.txt")
    return {"SRV": srv, "SRV+Double": mult2, "SRV+Double+Triple": mult23}


def parse_one_file(path: Path, run_type: str):
    txt = path.read_text(encoding="utf-8", errors="ignore")
    rep_m = re.search(r"_(\d+)_", path.name); rep = int(rep_m.group(1)) if rep_m else np.nan

    # p-value
    p = re.search(r"Likelihood ratio test\s+\**p\s*=\s*([0-9.]+)\**", txt); pval = float(p.group(1)) if p else np.nan

    # Block spans (for context mapping)
    alt_m = re.search(r"(Fitting the alternative model.*?)(?=\n###|\n----|\Z)", txt, re.S)
    null_m = re.search(r"(Fitting the null.*?)(?=\n###|\n----|\Z)", txt, re.S)
    alt_block = alt_m.group(1) if alt_m else ""
    null_block = null_m.group(1) if null_m else ""
    alt_span = (alt_m.start(1), alt_m.end(1)) if alt_m else None
    null_span = (null_m.start(1), null_m.end(1)) if null_m else None

    # K & metrics
    mk = re.search(r"Relaxation/intensification parameter\s*\(K\)\s*=\s*([0-9.]+)", alt_block); K = float(mk.group(1)) if mk else np.nan
    m_alt = re.search(r"\* Log\(L\)\s*=\s*([-\d.]+),\s*AIC-c\s*=\s*([-\d.]+)\s*\((\d+)\s+estimated parameters\)", alt_block)
    alt_logL, alt_AICc, alt_k = (float(m_alt.group(1)), float(m_alt.group(2)), int(m_alt.group(3))) if m_alt else (np.nan,np.nan,np.nan)
    m_null = re.search(r"\* Log\(L\)\s*=\s*([-\d.]+),\s*AIC-c\s*=\s*([-\d.]+)\s*\((\d+)\s+estimated parameters\)", null_block)
    null_logL, null_AICc, null_k = (float(m_null.group(1)), float(m_null.group(2)), int(m_null.group(3))) if m_null else (np.nan,np.nan,np.nan)

    # Global dN/dS
    mglob = re.search(r"(?:full codon model|Improving branch lengths).*?Reference\*\s*=\s*([0-9.]+).*?Test\*\s*=\s*([0-9.]+)", txt, re.S)
    glob_ref = float(mglob.group(1)) if mglob else np.nan
    glob_test = float(mglob.group(2)) if mglob else np.nan

    # Warnings
    warn_bits = []
    if re.search(r"Potential convergence issues|flat likelihood", txt): warn_bits.append("Convergence/flat likelihood")
    if re.search(r"Collapsed rate class", txt): warn_bits.append("Collapsed rate class")
    warnings = "; ".join(warn_bits)

    # Omega classes
    omega_rows = []
    def parse_omega(block, model_label):
        for set_label in ("test", "reference"):
            mtab = re.search(rf"^\*\s*The following rate distribution.*?\*\*{set_label}\*\* branches\s*\n((?:\|.*\n)+)", block, re.M|re.S)
            if mtab:
                for line in mtab.group(1).strip().splitlines():
                    cells = [c.strip() for c in line.strip().strip('|').split('|')]
                    if len(cells) >= 3:
                        try:
                            omega = float(cells[1]); prop = float(cells[2])
                            omega_rows.append({"model": model_label, "branch_set": set_label, "mode": cells[0], "omega": omega, "proportion_percent": prop})
                        except ValueError:
                            pass
        mtab2 = re.search(r"^\*\s*The following rate distribution for test/reference branches was inferred\s*\n((?:\|.*\n)+)", block, re.M|re.S)
        if mtab2:
            for line in mtab2.group(1).strip().splitlines():
                cells = [c.strip() for c in line.strip().strip('|').split('|')]
                if len(cells) >= 3:
                    try:
                        omega = float(cells[1]); prop = float(cells[2])
                        for set_label in ("test","reference"):
                            omega_rows.append({"model": model_label, "branch_set": set_label, "mode": cells[0], "omega": omega, "proportion_percent": prop})
                    except ValueError:
                        pass
    parse_omega(alt_block, "alternative")
    parse_omega(null_block, "null")

    # SRV classes with class_id
    srv_rows = []
    def parse_srv(block, model_label):
        m = re.search(r"^\*\s*The following rate distribution for site-to-site \*\*synonymous\*\* rate variation was inferred\s*\n((?:\|.*\n)+)", block, re.M|re.S)
        if not m: return
        class_id = 1
        for line in m.group(1).strip().splitlines():
            cells = [c.strip() for c in line.strip().strip('|').split('|')]
            if len(cells) >= 2:
                try:
                    rate = float(cells[0]); prop = float(cells[1])
                    srv_rows.append({"model": model_label, "class_id": class_id, "syn_rate": rate, "proportion_percent": prop})
                    class_id += 1
                except ValueError:
                    pass
    parse_srv(alt_block, "alternative")
    parse_srv(null_block, "null")

    # Multihit rates with span-based model inference
    mh_rows = []
    def parse_mh_in(block, model_label):
        for set_label in ("reference","test"):
            for hit, pat in (("2H", r"\*\s*2H rate for \*\*{set}\*\*:\s*([0-9.]+)"),
                             ("3H", r"\*\s*3H rate for \*\*{set}\*\*:\s*([0-9.]+)")):
                mm = re.search(pat.format(set=set_label), block)
                if mm:
                    mh_rows.append({"model": model_label, "branch_set": set_label, "hit_class": hit, "rate": float(mm.group(1))})
    parse_mh_in(alt_block, "alternative")
    parse_mh_in(null_block, "null")

    for mline in re.finditer(r"\*\s*([23])H rate for \*\*(reference|test)\*\*:\s*([0-9.]+)", txt):
        hit = f"{mline.group(1)}H"; bset = mline.group(2); rate = float(mline.group(3))
        pos = mline.start()
        # assign based on enclosing or nearest block
        model = "unspecified"
        if alt_span and alt_span[0] <= pos <= alt_span[1]:
            model = "alternative"
        elif null_span and null_span[0] <= pos <= null_span[1]:
            model = "null"
        else:
            if alt_span and null_span:
                model = "alternative" if pos < null_span[0] else "null"
            elif alt_span:
                model = "alternative"
            elif null_span:
                model = "null"
        if not any((r["hit_class"]==hit and r["branch_set"]==bset and abs(r["rate"]-rate)<1e-12 and r["model"]==model) for r in mh_rows):
            mh_rows.append({"model": model, "branch_set": bset, "hit_class": hit, "rate": rate})

    summary_row = {"replicate": rep, "run_type": run_type, "p_value": pval,
                   "K_alt": K, "logL_alt": alt_logL, "AICc_alt": alt_AICc, "params_alt": alt_k,
                   "logL_null": null_logL, "AICc_null": null_AICc, "params_null": null_k,
                   "global_dNdS_reference": glob_ref, "global_dNdS_test": glob_test, "warnings": warnings}
    return summary_row, omega_rows, srv_rows, mh_rows
def parse_all(base: Path, prefix: str):
    files_by_type = collect_files(base, prefix)
    summary, omega_all, srv_all, mh_all = [], [], [], []
    for run_type, files in files_by_type.items():
        for f in files:
            s, o, srv, mh = parse_one_file(f, run_type)
            for r in o: r.update({"run_type": run_type, "replicate": s["replicate"]})
            for r in srv: r.update({"run_type": run_type, "replicate": s["replicate"]})
            for r in mh: r.update({"run_type": run_type, "replicate": s["replicate"]})
            summary.append(s); omega_all += o; srv_all += srv; mh_all += mh

    def _sort(df, cols, empty_cols):
        if not len(df): return pd.DataFrame(columns=empty_cols)
        return pd.DataFrame(df).sort_values(cols).reset_index(drop=True)

    summary_df = pd.DataFrame(summary).sort_values(["run_type","replicate"]).reset_index(drop=True)

    omega_df = _sort(omega_all, ["run_type","replicate","model","branch_set","mode"],
                     ["run_type","replicate","model","branch_set","mode","omega","proportion_percent"])
    srv_df = _sort(srv_all, ["run_type","replicate","model","class_id"],
                   ["run_type","replicate","model","class_id","syn_rate","proportion_percent"])
    mh_df = _sort(mh_all, ["run_type","replicate","model","branch_set","hit_class"],
                  ["run_type","replicate","model","branch_set","hit_class","rate"])
    if not mh_df.empty:
        mh_df['model'] = mh_df['model'].fillna('unspecified').astype(str)

    return summary_df, omega_df, srv_df, mh_df

def stats_by_run_type(summary_df: pd.DataFrame, mh_df: pd.DataFrame):
    def medminmax(s):
        s = pd.to_numeric(s, errors="coerce").dropna()
        if s.empty: return (np.nan, np.nan, np.nan)
        return (float(s.median()), float(s.min()), float(s.max()))
    rows = []
    for rt, g in summary_df.groupby("run_type"):
        p_med, p_min, p_max = medminmax(g["p_value"])
        k_med, k_min, k_max = medminmax(g["K_alt"])
        wr_med, wr_min, wr_max = medminmax(g["global_dNdS_reference"])
        wt_med, wt_min, wt_max = medminmax(g["global_dNdS_test"])
        warn_n = int((g["warnings"].fillna("").astype(str)!="").sum())
        row = {"run_type": rt, "replicates": int(g["replicate"].nunique()),
               "p_median": p_med, "p_min": p_min, "p_max": p_max,
               "K_median": k_med, "K_min": k_min, "K_max": k_max,
               "omega_ref_median": wr_med, "omega_test_median": wt_med,
               "warnings_count": warn_n}
        sub = mh_df[mh_df["run_type"]==rt]
        for hit in ("2H","3H"):
            for bset in ("reference","test"):
                s = pd.to_numeric(sub[(sub["hit_class"]==hit)&(sub["branch_set"]==bset)]["rate"], errors="coerce").dropna()
                row[f"{hit}_{bset}_median"] = float(s.median()) if not s.empty else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values("run_type")

def make_pdf(out_pdf: Path, summary_df: pd.DataFrame, stats_df: pd.DataFrame,
             mh_df: pd.DataFrame, omega_df: pd.DataFrame, srv_df: pd.DataFrame):
    with PdfPages(out_pdf) as pdf:
        page_no = 1
        # Cover
        fig = plt.figure(figsize=(PAGE_W, PAGE_H)); _adjust_layout(fig)
        ax = fig.add_subplot(111); ax.axis('off')
        ax.text(0.5, 0.84, "RELAX Combined Report", ha='center', va='center', fontsize=18)
        ax.text(0.5, 0.79, "SRV vs SRV+Double vs SRV+Double+Triple", ha='center', va='center', fontsize=12)
        ax.text(0.1, 0.72, "Includes:", fontsize=10, ha='left')
        ax.text(0.12, 0.68, "• Run-type summary", fontsize=9, ha='left')
        ax.text(0.12, 0.64, "• Boxplots: K, p-values, 2H/3H", fontsize=9, ha='left')
        ax.text(0.12, 0.60, "• Full omega/SRV/multihit tables (column-paged)", fontsize=9, ha='left')
        _page_number(fig, page_no); pdf.savefig(fig); plt.close(fig); page_no += 1

        # Run-type summary table
        fig, ax = plt.subplots(figsize=(PAGE_W, PAGE_H)); _adjust_layout(fig, ax)
        cols = ["run_type","replicates","p_median","p_min","p_max","K_median","K_min","K_max",
                "omega_ref_median","omega_test_median","warnings_count","2H_reference_median","2H_test_median","3H_reference_median","3H_test_median"]
        for c in cols:
            if c not in stats_df.columns: stats_df[c] = np.nan
        draw_table(ax, stats_df[cols].round(4), "Run-type summary statistics", col_wrap={"run_type":10}, fontsize=8, scale=(1.55,1.55), header_wrap=14)
        _page_number(fig, page_no); pdf.savefig(fig); plt.close(fig); page_no += 1

        # Boxplots K and p
        order = list(stats_df["run_type"])
        fig, ax = plt.subplots(figsize=(PAGE_W, PAGE_H*0.66)); _adjust_layout(fig, ax)
        data = [pd.to_numeric(summary_df[summary_df["run_type"]==rt]["K_alt"], errors="coerce").dropna().values for rt in order]
        ax.boxplot(data, labels=order); ax.set_title("K (alternative) by run type"); ax.set_ylabel("K_alt")
        _page_number(fig, page_no); pdf.savefig(fig); plt.close(fig); page_no += 1

        fig, ax = plt.subplots(figsize=(PAGE_W, PAGE_H*0.66)); _adjust_layout(fig, ax)
        data = [pd.to_numeric(summary_df[summary_df["run_type"]==rt]["p_value"], errors="coerce").dropna().values for rt in order]
        ax.boxplot(data, labels=order); ax.set_title("RELAX LRT p-values by run type"); ax.set_ylabel("p-value"); ax.set_ylim(0,1)
        _page_number(fig, page_no); pdf.savefig(fig); plt.close(fig); page_no += 1

        # Boxplots for 2H/3H
        for hit in ("2H","3H"):
            mh_sub = mh_df[mh_df["hit_class"]==hit]
            if mh_sub.empty:
                fig, ax = plt.subplots(figsize=(PAGE_W, PAGE_H)); _adjust_layout(fig, ax); ax.axis('off')
                ax.text(0.5, 0.5, f"No {hit} rates available", ha='center', va='center')
                _page_number(fig, page_no); pdf.savefig(fig); plt.close(fig); page_no += 1
            else:
                for bset in ("reference","test"):
                    fig, ax = plt.subplots(figsize=(PAGE_W, PAGE_H*0.66)); _adjust_layout(fig, ax)
                    data = [pd.to_numeric(mh_sub[(mh_sub['run_type']==rt)&(mh_sub['branch_set']==bset)]['rate'], errors='coerce').dropna().values for rt in order]
                    clean = [d if len(d)>0 else np.array([np.nan]) for d in data]
                    ax.boxplot(clean, labels=order); ax.set_title(f"{hit} rates ({bset}) by run type"); ax.set_ylabel(f"{hit} rate")
                    _page_number(fig, page_no); pdf.savefig(fig); plt.close(fig); page_no += 1

        # Replicate-level summary by run_type
        cols_r = ["run_type","replicate","p_value","K_alt","global_dNdS_reference","global_dNdS_test","warnings"]
        for rt in order:
            sub = summary_df[summary_df["run_type"]==rt][cols_r].round(6)
            rows_per_page = 26
            pages = int(math.ceil(len(sub)/rows_per_page)) if len(sub)>0 else 1
            for i in range(pages):
                fig, ax = plt.subplots(figsize=(PAGE_W, PAGE_H)); _adjust_layout(fig, ax)
                chunk = sub.iloc[i*rows_per_page:(i+1)*rows_per_page]
                draw_table(ax, chunk, f"Replicate-level summary — {rt} (page {i+1}/{pages})",
                           col_wrap={"run_type":10, "warnings":28}, fontsize=7, scale=(1.5,1.5), header_wrap=14)
                _page_number(fig, page_no); pdf.savefig(fig); plt.close(fig); page_no += 1

        # Full omega, SRV, multihit tables with column paging
        page_ref = [page_no]
        render_wide_table(pdf, "Ω-class distributions (all)", omega_df, cols_per_page=15, col_wrap={"mode":10}, fontsize=7, scale=(1.45,1.45), page_no_ref=page_ref)
        page_no = page_ref[0]
        page_ref = [page_no]
        render_wide_table(pdf, "SRV rate classes (all)", srv_df, cols_per_page=15, col_wrap={"run_type":10}, fontsize=7, scale=(1.45,1.45), page_no_ref=page_ref)
        page_no = page_ref[0]
        page_ref = [page_no]
        render_wide_table(pdf, "Multihit (2H/3H) rates (all)", mh_df, cols_per_page=15, col_wrap={"branch_set":10}, fontsize=7, scale=(1.45,1.45), page_no_ref=page_ref)
        page_no = page_ref[0]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=".", help="Directory with HyPhy outputs")
    ap.add_argument("--prefix", required=True, help="File prefix, e.g., Camelus_bactrianus")
    args = ap.parse_args()

    base = Path(args.base)
    summary_df, omega_df, srv_df, mh_df = parse_all(base, args.prefix)

    # Save CSVs
    (base / f"{args.prefix}_RELAX_ALL_summary.csv").write_text(summary_df.to_csv(index=False), encoding="utf-8")
    (base / f"{args.prefix}_RELAX_ALL_omega_classes.csv").write_text(omega_df.to_csv(index=False), encoding="utf-8")
    (base / f"{args.prefix}_RELAX_ALL_syn_rate_classes.csv").write_text(srv_df.to_csv(index=False), encoding="utf-8")
    (base / f"{args.prefix}_RELAX_ALL_multihit_rates.csv").write_text(mh_df.to_csv(index=False), encoding="utf-8")

    # Extra: combined classes (omega + SRV) long table for convenience
    comb_rows = []
    for _, r in omega_df.iterrows():
        comb_rows.append({"run_type": r.get("run_type"), "replicate": r.get("replicate"), "model": r.get("model"),
                           "branch_set": r.get("branch_set"), "class_type": "omega", "class_id": r.get("mode"),
                           "omega": r.get("omega"), "syn_rate": None, "proportion_percent": r.get("proportion_percent")})
    for _, r in srv_df.iterrows():
        comb_rows.append({"run_type": r.get("run_type"), "replicate": r.get("replicate"), "model": r.get("model"),
                           "branch_set": None, "class_type": "srv", "class_id": r.get("class_id"),
                           "omega": None, "syn_rate": r.get("syn_rate"), "proportion_percent": r.get("proportion_percent")})
    comb_df = pd.DataFrame(comb_rows)
    (base / f"{args.prefix}_RELAX_ALL_classes_long.csv").write_text(comb_df.to_csv(index=False), encoding="utf-8")

    # Stats and PDF
    stats_df = stats_by_run_type(summary_df, mh_df)
    make_pdf(base / f"{args.prefix}_RELAX_ALL_report.pdf", summary_df, stats_df, mh_df, omega_df, srv_df)

if __name__ == "__main__":
    main()
