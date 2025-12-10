#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from pathlib import Path
import itertools

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Ejecuta tests no paramétricos (Friedman + Wilcoxon pareado con FDR BH) "
            "sobre un fichero de métricas por muestra en formato ancho.\n\n"
            "Entrada esperada: TSV con columnas\n"
            "  - herramienta\n"
            "  - muestra\n"
            "  - una columna por métrica (Aitchison, JSD, BrayCurtis, Spearman, F1, ...)\n"
        )
    )
    ap.add_argument(
        "metricas_por_muestra",
        help="TSV con las métricas por muestra y herramienta (formato ancho).",
    )
    ap.add_argument(
        "--outdir",
        default="tests_metricas",
        help="Directorio de salida (default: tests_metricas).",
    )
    return ap.parse_args()


def fdr_bh(pvals, alpha=0.05):
    """
    Benjamini-Hochberg FDR.

    pvals: array-like de p-values (sin NaNs).
    Devuelve:
      - reject: array bool (True si significativo con FDR<=alpha)
      - p_adj: p-values ajustados.
    """
    p = np.asarray(pvals, dtype=float)
    m = p.size
    if m == 0:
        return np.array([], dtype=bool), np.array([], dtype=float)

    # Ordenar p-values
    order = np.argsort(p)
    ranked_p = p[order]

    # p ajustados según BH
    adj = np.empty(m, dtype=float)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        rank = i + 1
        val = ranked_p[i] * m / rank
        if val < prev:
            prev = val
        adj[i] = prev

    # Volver al orden original
    p_adj = np.empty(m, dtype=float)
    p_adj[order] = adj

    reject = p_adj <= alpha
    return reject, p_adj


def main():
    args = parse_args()
    in_path = Path(args.metricas_por_muestra)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Leer tabla de métricas (ancha)
    try:
        df = pd.read_csv(in_path, sep="\t")
    except Exception as e:
        print(f"[ERROR] No se pudo leer {in_path}: {e}", file=sys.stderr)
        sys.exit(1)

    required = {"herramienta", "muestra"}
    if not required.issubset(df.columns):
        print(
            f"[ERROR] Faltan columnas obligatorias. Se esperaban al menos: {required}. "
            f"Columnas encontradas: {list(df.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Columnas de métricas = todo lo que no sea herramienta ni muestra
    metric_cols = [c for c in df.columns if c not in ("herramienta", "muestra")]
    if not metric_cols:
        print(
            "[ERROR] No se detectaron columnas de métricas aparte de 'herramienta' y 'muestra'.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"[INFO] Métricas detectadas: {metric_cols}", file=sys.stderr)

    # ---------- TEST DE FRIEDMAN POR MÉTRICA ----------
    friedman_rows = []

    for met in metric_cols:
        sub = df[["muestra", "herramienta", met]].copy()

        # Tabla muestras × herramientas
        wide = sub.pivot(index="muestra", columns="herramienta", values=met)

        # Eliminar muestras con NaN en alguna herramienta
        wide = wide.dropna(axis=0, how="any")

        # Necesitamos al menos 2 herramientas y varias muestras
        if wide.shape[0] < 2 or wide.shape[1] < 2:
            print(
                f"[WARN] Métrica {met}: no hay datos suficientes para Friedman "
                f"(muestras={wide.shape[0]}, herramientas={wide.shape[1]}).",
                file=sys.stderr,
            )
            continue

        try:
            stats = [wide[col].values for col in wide.columns]
            stat, p = friedmanchisquare(*stats)
        except Exception as e:
            print(f"[WARN] Error en Friedman para {met}: {e}", file=sys.stderr)
            continue

        friedman_rows.append(
            {"metrica": met, "friedman_stat": stat, "friedman_p": p}
        )

    friedman_df = pd.DataFrame(friedman_rows)
    friedman_out = outdir / "friedman_metricas.tsv"
    friedman_df.to_csv(friedman_out, sep="\t", index=False)
    print(f"[OK] Test de Friedman guardado en: {friedman_out}", file=sys.stderr)

    # ---------- WILCOXON PAREADO + FDR (BH) POR MÉTRICA ----------
    wilcoxon_all = []

    for met in metric_cols:
        sub = df[["muestra", "herramienta", met]].copy()
        wide = sub.pivot(index="muestra", columns="herramienta", values=met)
        wide = wide.dropna(axis=0, how="any")

        tools = list(wide.columns)
        if len(tools) < 2:
            print(
                f"[WARN] Métrica {met}: menos de 2 herramientas tras filtrar, "
                f"se omite Wilcoxon.",
                file=sys.stderr,
            )
            continue

        pair_rows = []
        p_values = []

        for h1, h2 in itertools.combinations(tools, 2):
            v1 = wide[h1].values
            v2 = wide[h2].values

            if len(v1) < 2:
                continue

            try:
                stat, p = wilcoxon(
                    v1, v2, zero_method="wilcox", alternative="two-sided"
                )
            except ValueError as e:
                # Suele pasar si todas las diferencias son cero
                print(
                    f"[WARN] Wilcoxon falló para {met} ({h1} vs {h2}): {e}",
                    file=sys.stderr,
                )
                stat, p = float("nan"), float("nan")

            pair_rows.append(
                {
                    "metrica": met,
                    "h1": h1,
                    "h2": h2,
                    "stat": stat,
                    "p": p,
                }
            )
            p_values.append(p)

        if not pair_rows:
            continue

        # FDR Benjamini–Hochberg por métrica
        p_series = pd.Series([r["p"] for r in pair_rows], dtype="float64")
        valid_mask = p_series.notna()

        if valid_mask.any():
            p_valid = p_series[valid_mask].values
            reject, p_adj = fdr_bh(p_valid, alpha=0.05)

            idx_valid = [i for i, ok in enumerate(valid_mask) if ok]
            for idx_local, (rej, padj) in enumerate(zip(reject, p_adj)):
                global_idx = idx_valid[idx_local]
                pair_rows[global_idx]["p_fdr"] = padj
                pair_rows[global_idx]["significativo"] = bool(rej)

        # Rellenar NaNs en los que no se ajustaron
        for r in pair_rows:
            if "p_fdr" not in r:
                r["p_fdr"] = float("nan")
                r["significativo"] = False

        wilcoxon_all.extend(pair_rows)

    wilcoxon_df = pd.DataFrame(wilcoxon_all)
    wilcoxon_out = outdir / "wilcoxon_pares.tsv"
    wilcoxon_df.to_csv(wilcoxon_out, sep="\t", index=False)
    print(f"[OK] Tests de Wilcoxon pareado + FDR (BH) guardados en: {wilcoxon_out}", file=sys.stderr)


if __name__ == "__main__":
    main()

