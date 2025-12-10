#!/usr/bin/env python3
# Compara tu unified_abundance.clean.csv con el Excel del paper (Relative abundances (long))
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path
import sys

CLEAN = Path("results/taxpasta/unified_abundance.clean.csv")
# ⚠️ Ajusta la ruta al Excel en TU HPC (colócalo donde quieras y cambia aquí):
PAPER_XLSX = Path("../../RAW/ERP015409/sup_data_3_ERP015409.xlsx")

OUTD  = Path("results/metrics")
OUTD.mkdir(parents=True, exist_ok=True)

def bray_curtis(a, b):
    den = (a + b).sum()
    return np.nan if den <= 0 else np.abs(a - b).sum() / den

def main():
    if not CLEAN.exists():
        sys.exit(f"No existe {CLEAN}. Ejecuta primero el clean_unified.py.")
    if not PAPER_XLSX.exists():
        sys.exit(f"No encuentro el Excel del paper en {PAPER_XLSX}. Copia el .xlsx ahí o cambia la ruta en el script.")

    # 1) Tus datos (ya limpios a nivel genus)
    your = pd.read_csv(CLEAN, dtype={"sample_id": str, "profiler": str, "tax_name": str})
    your = your[["sample_id","profiler","tax_name","rel_abundance"]].copy()
    your["tax_lc"] = your["tax_name"].str.lower()

    # 2) Paper (hoja "Relative abundances (long)")
    try:
        ref = pd.read_excel(PAPER_XLSX, sheet_name="Relative abundances (long)")
    except ImportError as e:
        sys.exit("Pandas necesita 'openpyxl' para leer .xlsx. Instálalo con:\n"
                 "  python -m pip install --user openpyxl\n"
                 f"Error original: {e}")
    ref.columns = [c.strip().lower().replace(" ","_") for c in ref.columns]
    if "sample" in ref.columns: ref = ref.rename(columns={"sample":"sample_id"})
    if "genus"  in ref.columns: ref = ref.rename(columns={"genus":"tax_name"})
    if "relative_abundance" not in ref.columns and "abundance" in ref.columns:
        ref = ref.rename(columns={"abundance":"relative_abundance"})
    if ref["relative_abundance"].max() > 1.0:
        ref["relative_abundance"] = ref["relative_abundance"] / 100.0
    ref = ref[["sample_id","tax_name","relative_abundance"]].copy()
    ref["tax_lc"] = ref["tax_name"].astype(str).str.lower()

    # 3) Emparejar por (sample_id, taxón)
    merged = (your.groupby(["sample_id","profiler","tax_lc"], as_index=False)["rel_abundance"].sum()
                   .merge(ref.groupby(["sample_id","tax_lc"], as_index=False)["relative_abundance"].sum(),
                          on=["sample_id","tax_lc"], how="inner"))

    print(f"Muestras en común: {merged['sample_id'].nunique()} | pares taxón-muestra: {len(merged)}")

    # 4) Métricas globales por perfilador
    rows = []
    for p in sorted(merged["profiler"].unique()):
        sub = merged[merged["profiler"]==p]
        rho = spearmanr(sub["rel_abundance"], sub["relative_abundance"]).correlation
        bc  = bray_curtis(sub["rel_abundance"].values, sub["relative_abundance"].values)
        rows.append({"profiler": p, "spearman_rho": rho, "bray_curtis": bc, "n_pairs": len(sub)})
    scores = pd.DataFrame(rows).sort_values("spearman_rho", ascending=False)
    scores.to_csv(OUTD/"global_scores_by_profiler.csv", index=False)
    print("\n== Global scores by profiler ==")
    print(scores.to_string(index=False))

    # 5) Top-10 overlap por muestra y perfilador
    def topk_overlap(df_y, df_r, k=10):
        s1 = set(df_y.sort_values("rel_abundance", ascending=False).head(k)["tax_lc"])
        s2 = set(df_r.sort_values("relative_abundance", ascending=False).head(k)["tax_lc"])
        return np.nan if (len(s1)==0 and len(s2)==0) else len(s1 & s2) / max(1, len(s1 | s2))

    topk_rows = []
    for p in sorted(merged["profiler"].unique()):
        for s in merged["sample_id"].unique():
            y = merged[(merged["profiler"]==p) & (merged["sample_id"]==s)][["tax_lc","rel_abundance"]]
            r = merged[merged["sample_id"]==s][["tax_lc","relative_abundance"]].drop_duplicates()
            if len(y)==0 or len(r)==0:
                continue
            topk_rows.append({"profiler": p, "sample_id": s, "top10_overlap": topk_overlap(y, r, 10)})
    topk = pd.DataFrame(topk_rows)
    topk.to_csv(OUTD/"top10_overlap_by_sample_profiler.csv", index=False)
    print("\n== Top-10 overlap (primeras filas) ==")
    if not topk.empty:
        print(topk.head().to_string(index=False))
    else:
        print("Sin filas (revisa emparejamiento de sample_id/tax_name).")

if __name__ == "__main__":
    main()
