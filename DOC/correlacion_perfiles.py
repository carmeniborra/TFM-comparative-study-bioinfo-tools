#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from scipy.stats import pearsonr, spearmanr, friedmanchisquare, wilcoxon

# ---------- CARGA Y ALINEACIÓN ----------

def cargar_tabla(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df.iloc[:, 1:]  # quitar columna de abundancia total


def alinear(ref, tool):
    genera = ref.index.intersection(tool.index)
    samples = ref.columns.intersection(tool.columns)
    ref = ref.loc[genera, samples].sort_index().sort_index(axis=1)
    tool = tool.loc[genera, samples].sort_index().sort_index(axis=1)
    return ref, tool


# ---------- BENJAMINI–HOCHBERG ----------

def fdr_bh(pvals, alpha=0.05):
    pvals = np.array(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(1, n+1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj[adj>1] = 1.0
    p_fdr = np.empty_like(adj)
    p_fdr[order] = adj
    return p_fdr, p_fdr < alpha


# ---------- MAIN ----------

def main():
    parser = argparse.ArgumentParser(
        description="Correlación Pearson y Spearman entre GS y herramientas."
    )
    parser.add_argument("ref_tsv", help="TSV GS")
    parser.add_argument("carpeta_tools", help="Carpeta con TSVs herramienta")
    args = parser.parse_args()

    ref_path = Path(args.ref_tsv)
    tools_folder = Path(args.carpeta_tools)

    ref = cargar_tabla(ref_path)

    outdir = Path.cwd() / "correlacion_perfiles"
    outdir.mkdir(exist_ok=True)

    filas = []

    for tool_path in sorted(tools_folder.glob("*.tsv")):

        if tool_path.resolve() == ref_path.resolve():
            continue

        tool = cargar_tabla(tool_path)
        ref_al, tool_al = alinear(ref, tool)

        for sample in ref_al.columns:
            x = ref_al[sample].values
            y = tool_al[sample].values

            # Pearson
            r_pearson, p1 = pearsonr(x, y)

            # Spearman
            r_spearman, p2 = spearmanr(x, y)

            filas.append({
                "herramienta": tool_path.name,
                "muestra": sample,
                "pearson_r": r_pearson,
                "spearman_r": r_spearman
            })

    df = pd.DataFrame(filas)
    df.to_csv(outdir / "correlaciones_por_muestra.tsv", sep="\t", index=False)

    # ---------- RESUMEN POR HERRAMIENTA ----------
    resumen = (
        df.groupby("herramienta")[["pearson_r", "spearman_r"]]
          .agg(["mean", "median", "min", "max"])
    )
    resumen.to_csv(outdir / "resumen_correlacion.tsv", sep="\t")

    # ---------- FRIEDMAN ----------
    def friedman(metric):
        tabla = df.pivot(index="muestra", columns="herramienta", values=metric).dropna()
        datos = [tabla[h].values for h in tabla.columns]
        stat, p = friedmanchisquare(*datos)
        return stat, p, tabla

    stat_p, p_p, tabpear = friedman("pearson_r")
    stat_s, p_s, tabspear = friedman("spearman_r")

    pd.DataFrame([
        {"metrica":"pearson_r", "friedman_stat": stat_p, "p": p_p},
        {"metrica":"spearman_r", "friedman_stat": stat_s, "p": p_s}
    ]).to_csv(outdir / "friedman_correlacion.tsv", sep="\t", index=False)

    # ---------- WILCOXON + FDR ----------
    pairs = []
    for metric, tabla in zip(["pearson_r","spearman_r"], [tabpear, tabspear]):
        herramientas = tabla.columns
        pvals = []
        comp_tmp = []

        for i in range(len(herramientas)):
            for j in range(i+1, len(herramientas)):
                h1, h2 = herramientas[i], herramientas[j]
                x, y = tabla[h1].values, tabla[h2].values

                try:
                    stat, p = wilcoxon(x, y)
                except:
                    stat, p = np.nan, np.nan

                comp_tmp.append({
                    "metrica": metric,
                    "h1": h1,
                    "h2": h2,
                    "stat": stat,
                    "p": p
                })
                pvals.append(p if not np.isnan(p) else 1.0)

        pfdr, sig = fdr_bh(pvals)
        for k,r in enumerate(comp_tmp):
            r["p_fdr"] = pfdr[k]
            r["significativo"] = sig[k]
            pairs.append(r)

    pd.DataFrame(pairs).to_csv(outdir / "wilcoxon_correlacion.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

