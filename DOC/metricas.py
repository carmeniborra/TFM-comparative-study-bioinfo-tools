#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis
from scipy.stats import spearmanr

# --------------------
# Funciones auxiliares
# --------------------

def load_abundance_table(path):
    """
    Lee una tabla de abundancias con:
    - fila = taxón (género)
    - col0 = algo tipo "sum_abundance"
    - resto de columnas = muestras
    Devuelve DataFrame con:
    - index = taxón
    - columns = solo muestras
    """
    df = pd.read_csv(path, sep="\t", header=0, index_col=0)
    # asumimos que la primera columna es sum_abundance y la quitamos
    if df.shape[1] > 1:
        df = df.iloc[:, 1:]
    df = df.astype(float).fillna(0.0)
    return df


def closure(x):
    """Cierra un vector composicional para que sume 1."""
    x = np.asarray(x, dtype=float)
    s = x.sum()
    if s == 0:
        return np.ones_like(x) / len(x)
    return x / s


def clr(x, pseudocount=1e-6):
    """
    Transformación CLR con pseudoconteo.
    x: vector composicional (no se asume cerrado).
    """
    x = np.asarray(x, dtype=float)
    x = x + pseudocount
    x = closure(x)
    g_mean = np.exp(np.mean(np.log(x)))
    return np.log(x / g_mean)


def aitchison_distance(p, q, pseudocount=1e-6):
    """Distancia de Aitchison entre dos composiciones p y q."""
    clr_p = clr(p, pseudocount=pseudocount)
    clr_q = clr(q, pseudocount=pseudocount)
    return np.linalg.norm(clr_p - clr_q)


def kl_divergence(p, q):
    """KL(p || q) con protección a ceros. Se asume p y q cerrados."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    return np.sum(p * np.log(p / q))


def jensen_shannon_divergence(p, q, pseudocount=1e-6):
    """
    Jensen–Shannon Divergence entre dos composiciones p y q.
    Devuelve JSD en unidades de log natural.
    """
    p = np.asarray(p, dtype=float) + pseudocount
    q = np.asarray(q, dtype=float) + pseudocount
    p = closure(p)
    q = closure(q)
    m = 0.5 * (p + q)
    jsd = 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)
    return jsd


def f1_score_binary(y_true, y_pred):
    """
    F1-score de presencia/ausencia:
    - presente si abundancia > 0
    """
    y_true = np.asarray(y_true, dtype=float) > 0
    y_pred = np.asarray(y_pred, dtype=float) > 0

    tp = np.logical_and(y_true, y_pred).sum()
    fp = np.logical_and(~y_true, y_pred).sum()
    fn = np.logical_and(y_true, ~y_pred).sum()

    denom = 2 * tp + fp + fn
    if denom == 0:
        return np.nan
    return 2 * tp / denom


def normalize_tool_name(filename):
    """
    Normaliza el nombre de la herramienta a partir del nombre de fichero.
    Ejemplos:
    - 'bracken_db2.tsv.genus.tsv.study.aligned.tsv' -> 'bracken_db2'
    - 'kaiju_db5.tsv.genus.tsv.study.aligned.tsv'  -> 'kaiju_db5'
    - 'gold_standard.aligned.tsv'                  -> 'gold_standard'
    """
    # nos quedamos con lo que hay antes del primer ".tsv"
    if ".tsv" in filename:
        return filename.split(".tsv")[0]
    else:
        # por si acaso, quitamos extensión simple .tsv / .txt
        return os.path.splitext(filename)[0]


# --------------------
# Cálculo de métricas
# --------------------

def compute_metrics_for_tool(gs_df, tool_df, herramienta_name):
    """
    Calcula todas las métricas por muestra para una herramienta.
    Devuelve lista de dicts.
    """
    # Intersección de géneros y muestras
    common_genera = gs_df.index.intersection(tool_df.index)
    common_samples = gs_df.columns.intersection(tool_df.columns)

    gs_sub = gs_df.loc[common_genera, common_samples]
    tool_sub = tool_df.loc[common_genera, common_samples]

    resultados = []

    for sample in common_samples:
        ref = gs_sub[sample].values.astype(float)
        est = tool_sub[sample].values.astype(float)

        # Cerramos por si acaso
        ref_c = closure(ref)
        est_c = closure(est)

        # Aitchison
        aitch = aitchison_distance(ref_c, est_c)

        # JSD
        jsd = jensen_shannon_divergence(ref_c, est_c)

        # Bray-Curtis
        bc = braycurtis(ref_c, est_c)

        # Spearman
        rho, _ = spearmanr(ref, est)

        # F1 (presencia/ausencia)
        f1 = f1_score_binary(ref, est)

        resultados.append({
            "herramienta": herramienta_name,
            "muestra": sample,
            "Aitchison": aitch,
            "JSD": jsd,
            "BrayCurtis": bc,
            "Spearman": rho,
            "F1": f1
        })

    return resultados


def main():
    parser = argparse.ArgumentParser(
        description="Calcular métricas por muestra (Aitchison, JSD, Bray-Curtis, Spearman, F1)."
    )
    parser.add_argument(
        "--gs",
        required=True,
        help="Tabla de referencia (GS) en formato TSV (por ejemplo gold_standard.aligned.tsv)."
    )
    parser.add_argument(
        "--tools_dir",
        required=True,
        help="Carpeta con las tablas de las herramientas (TSV). Puede contener también el GS."
    )
    parser.add_argument(
        "--out",
        default="metricas_por_muestra.tsv",
        help="Fichero de salida con todas las métricas por muestra."
    )

    args = parser.parse_args()

    # Cargamos GS
    gs_df = load_abundance_table(args.gs)
    gs_basename = os.path.basename(args.gs)

    all_results = []

    for fname in sorted(os.listdir(args.tools_dir)):
        if not fname.endswith(".tsv"):
            continue

        # Saltamos el GS si está en la misma carpeta
        if fname == gs_basename or "gold_standard" in fname:
            continue

        tool_path = os.path.join(args.tools_dir, fname)
        tool_name = normalize_tool_name(fname)

        print(f"Analizando {fname} -> herramienta '{tool_name}'...")
        tool_df = load_abundance_table(tool_path)

        res = compute_metrics_for_tool(gs_df, tool_df, tool_name)
        all_results.extend(res)

    if not all_results:
        print("No se han encontrado resultados. Revisa la carpeta de herramientas y el GS.")
        return

    df_out = pd.DataFrame(all_results)
    df_out.to_csv(args.out, sep="\t", index=False)
    print(f"✅ Métricas por muestra guardadas en {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
