#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
import argparse

from scipy.stats import friedmanchisquare, wilcoxon, shapiro

# ---------- MÉTRICAS ----------

def l1_distance(x, y):
    return np.sum(np.abs(x - y))

def bray_curtis(x, y):
    return l1_distance(x, y) / (np.sum(x + y))

def f1_score(x, y, threshold=0.001):
    ref_detect = x > threshold
    tool_detect = y > threshold
    
    TP = np.sum(ref_detect & tool_detect)
    FP = np.sum(~ref_detect & tool_detect)
    FN = np.sum(ref_detect & ~tool_detect)

    # Caso especial: nadie detecta nada
    if TP == 0 and FP == 0 and FN == 0:
        return 1.0

    if TP == 0 and (FP + FN) > 0:
        return 0.0
    
    return 2 * TP / (2 * TP + FP + FN)

def clr_transform(v, pseudocount=1e-6):
    v = v + pseudocount
    gm = np.exp(np.mean(np.log(v)))
    return np.log(v / gm)

def aitchison_distance(x, y):
    clr_x = clr_transform(x)
    clr_y = clr_transform(y)
    return np.linalg.norm(clr_x - clr_y)

def clr_matrix(df, pseudocount=1e-6):
    """
    Aplica CLR por columnas (muestras) a una matriz de composiciones:
    filas = géneros, columnas = muestras.
    """
    arr = df.values.astype(float) + pseudocount
    gm = np.exp(np.mean(np.log(arr), axis=0))    # media geométrica por columna
    clr = np.log(arr / gm)
    return pd.DataFrame(clr, index=df.index, columns=df.columns)


# ---------- CARGA Y ALINEACIÓN ----------

def cargar_tabla(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    # quitar columna de abundancia acumulada
    return df.iloc[:, 1:]

def alinear(ref, tool):
    genera = ref.index.intersection(tool.index)
    samples = ref.columns.intersection(tool.columns)
    ref = ref.loc[genera, samples].sort_index().sort_index(axis=1)
    tool = tool.loc[genera, samples].sort_index().sort_index(axis=1)
    return ref, tool


# ---------- FDR (Benjamini–Hochberg) ----------

def fdr_bh(pvals, alpha=0.05):
    """
    Implementación sencilla de Benjamini–Hochberg.
    Devuelve:
      reject: booleanos (significativo a FDR<alpha)
      p_fdr: p-valores ajustados
    """
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    if n == 0:
        return np.array([], dtype=bool), np.array([], dtype=float)

    order = np.argsort(pvals)
    ranked_p = pvals[order]

    adj = ranked_p * n / (np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj[adj > 1] = 1.0

    p_fdr = np.empty_like(adj)
    p_fdr[order] = adj

    reject = p_fdr < alpha
    return reject, p_fdr


# ---------- ESTADÍSTICA ENTRE HERRAMIENTAS ----------

def pruebas_entre_herramientas(df_metrics, output_folder):
    """
    df_metrics: DataFrame con columnas:
        herramienta, muestra, L1, BrayCurtis, F1, Aitchison
    output_folder: carpeta donde se escriben los resultados.
    """
    metrics = ["L1", "BrayCurtis", "F1", "Aitchison"]
    global_rows = []
    pairwise_rows = []

    for met in metrics:
        if met not in df_metrics.columns:
            continue

        tabla = df_metrics.pivot(index="muestra", columns="herramienta", values=met)
        tabla = tabla.dropna(axis=0, how="any")

        if tabla.shape[0] < 2 or tabla.shape[1] < 2:
            print(f"No hay suficientes datos para {met} (Friedman).")
            continue

        herramientas = list(tabla.columns)

        # ---- Friedman ----
        try:
            datos_friedman = [tabla[h].values for h in herramientas]
            fried_stat, fried_p = friedmanchisquare(*datos_friedman)
        except Exception as e:
            print(f"Problema con Friedman en {met}: {e}")
            fried_stat, fried_p = np.nan, np.nan

        global_rows.append({
            "metrica": met,
            "friedman_stat": fried_stat,
            "friedman_p": fried_p
        })

        # ---- Wilcoxon + FDR ----
        pares_tmp = []
        pvals = []

        for i in range(len(herramientas)):
            for j in range(i + 1, len(herramientas)):
                h1 = herramientas[i]
                h2 = herramientas[j]

                x = tabla[h1].values
                y = tabla[h2].values

                try:
                    stat, p = wilcoxon(x, y, alternative="two-sided")
                except Exception as e:
                    print(f"Problema con Wilcoxon ({met}, {h1} vs {h2}): {e}")
                    stat, p = np.nan, np.nan

                pares_tmp.append({
                    "metrica": met,
                    "herramienta_1": h1,
                    "herramienta_2": h2,
                    "stat": stat,
                    "p_val": p
                })
                pvals.append(p if not np.isnan(p) else 1.0)

        if pvals:
            reject, p_fdr = fdr_bh(pvals, alpha=0.05)
            for k, row in enumerate(pares_tmp):
                row["p_val_FDR"] = p_fdr[k]
                row["significativo_FDR_0.05"] = bool(reject[k])
                pairwise_rows.append(row)

    if global_rows:
        df_global = pd.DataFrame(global_rows)
        df_global.to_csv(output_folder / "global_tests_Friedman.tsv",
                         sep="\t", index=False)
        print(f"✅ Tests globales (Friedman) guardados en {output_folder / 'global_tests_Friedman.tsv'}")

    if pairwise_rows:
        df_pairs = pd.DataFrame(pairwise_rows)
        df_pairs.to_csv(output_folder / "pairwise_Wilcoxon_FDR.tsv",
                        sep="\t", index=False)
        print(f"✅ Comparaciones por pares guardadas en {output_folder / 'pairwise_Wilcoxon_FDR.tsv'}")


# ---------- PRUEBAS DE NORMALIDAD ----------

def pruebas_normalidad(df_metrics, output_folder):
    """
    Shapiro–Wilk por métrica y herramienta.
    """
    metrics = ["L1", "BrayCurtis", "F1", "Aitchison"]
    rows = []

    for met in metrics:
        if met not in df_metrics.columns:
            continue

        for herramienta in df_metrics["herramienta"].unique():
            vals = df_metrics.loc[df_metrics["herramienta"] == herramienta, met].dropna()

            if len(vals) < 3:
                stat, p = np.nan, np.nan
            else:
                try:
                    stat, p = shapiro(vals.values)
                except Exception as e:
                    print(f"Problema con Shapiro ({met}, {herramienta}): {e}")
                    stat, p = np.nan, np.nan

            rows.append({
                "metrica": met,
                "herramienta": herramienta,
                "n_muestras": len(vals),
                "shapiro_W": stat,
                "shapiro_p": p
            })

    if rows:
        df_norm = pd.DataFrame(rows)
        df_norm.to_csv(output_folder / "normality_tests_Shapiro.tsv",
                       sep="\t", index=False)
        print(f"✅ Pruebas de normalidad (Shapiro) guardadas en {output_folder / 'normality_tests_Shapiro.tsv'}")


# ---------- MAIN ----------

def main():
    parser = argparse.ArgumentParser(
        description=("Calcular distancias (L1, Bray-Curtis, F1, Aitchison) por muestra y herramienta, "
                     "comparar herramientas con Friedman y Wilcoxon + FDR, "
                     "evaluar normalidad con Shapiro–Wilk y guardar matrices CLR en carpetas separadas.")
    )
    parser.add_argument("ref_tsv", help="Tabla de referencia TSV")
    parser.add_argument("carpeta_tools", help="Carpeta con los perfiles TSV (entrada)")
    args = parser.parse_args()

    ref_path = Path(args.ref_tsv)
    tools_folder = Path(args.carpeta_tools)

    if not ref_path.is_file():
        print(f"ERROR: La tabla de referencia no existe: {ref_path}")
        return

    if not tools_folder.is_dir():
        print(f"ERROR: La carpeta de herramientas no existe: {tools_folder}")
        return

    # Carpeta raíz de salida: desde donde se ejecuta el script
    cwd = Path.cwd()
    output_root = cwd / "compare_metrics"
    output_root.mkdir(exist_ok=True)

    # Carpeta para matrices CLR al mismo nivel
    clr_dir = cwd / "CLR_matrices"
    clr_dir.mkdir(exist_ok=True)

    print(f"Salida principal en: {output_root}")
    print(f"Matrices CLR en:      {clr_dir}\n")

    # Cargar referencia
    ref = cargar_tabla(ref_path)

    # Buscar TSV de herramientas
    tsv_files = sorted(tools_folder.glob("*.tsv"))
    if not tsv_files:
        print(f"No se han encontrado .tsv en {tools_folder}")
        return

    resultados = []

    for tool_path in tsv_files:
        # Por si la referencia está en la misma carpeta que las herramientas
        if tool_path.resolve() == ref_path.resolve():
            continue

        print(f"Analizando {tool_path.name}...")
        tool = cargar_tabla(tool_path)
        ref_al, tool_al = alinear(ref, tool)

        # ---------- MATRICES CLR (alineadas) ----------
        clr_ref = clr_matrix(ref_al)
        clr_tool = clr_matrix(tool_al)

        clr_ref_file = clr_dir / f"{tool_path.stem}_CLR_reference.tsv"
        clr_tool_file = clr_dir / f"{tool_path.stem}_CLR_tool.tsv"

        clr_ref.to_csv(clr_ref_file, sep="\t")
        clr_tool.to_csv(clr_tool_file, sep="\t")

        # ---------- MÉTRICAS POR MUESTRA ----------
        for sample in ref_al.columns:
            x = ref_al[sample].values
            y = tool_al[sample].values

            L1 = l1_distance(x, y)
            BC = bray_curtis(x, y)
            F1 = f1_score(x, y)
            AD = aitchison_distance(x, y)

            resultados.append({
                "herramienta": tool_path.name,
                "muestra": sample,
                "L1": L1,
                "BrayCurtis": BC,
                "F1": F1,
                "Aitchison": AD
            })

    df_metrics = pd.DataFrame(resultados)

    # Guardar métricas por muestra
    out_metrics = output_root / "metricas_por_muestra.tsv"
    df_metrics.to_csv(out_metrics, sep="\t", index=False)
    print(f"\n✅ Métricas por muestra guardadas en {out_metrics}")

    # Pruebas de normalidad
    pruebas_normalidad(df_metrics, output_root)

    # Pruebas entre herramientas
    pruebas_entre_herramientas(df_metrics, output_root)


if __name__ == "__main__":
    main()

