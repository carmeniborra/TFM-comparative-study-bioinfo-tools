#!/usr/bin/env python3
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Para entornos sin display (HPC)
import matplotlib.pyplot as plt


def preparar_resumen(f_resumen):
    """
    Lee resumen_correlacion.tsv y devuelve un DataFrame con:
      - index = herramienta
      - columnas: pearson_mean, pearson_median, pearson_min, pearson_max,
                  spearman_mean, spearman_median, spearman_min, spearman_max
    Sea cual sea el formato original (doble cabecera, etc.).
    """
    df = pd.read_csv(f_resumen, sep="\t", header=0, index_col=0)

    # Si ya está en formato "bonito" y tiene pearson_median, devolvemos tal cual
    if "pearson_median" in df.columns:
        return df

    cols = df.columns.tolist()

    # NUEVA detección: columnas que empiezan por pearson_r / spearman_r
    hay_pearson = any(c.startswith("pearson_r") for c in cols)
    hay_spearman = any(c.startswith("spearman_r") for c in cols)

    if hay_pearson and hay_spearman:
        # Primera fila contiene los sufijos: mean, median, min, max...
        sufijos = df.iloc[0]

        nuevos_nombres = []
        for base, sufijo in zip(cols, sufijos):
            base_clean = base.strip()
            suf_clean = str(sufijo).strip()

            if base_clean.startswith("pearson"):
                pref = "pearson"
            elif base_clean.startswith("spearman"):
                pref = "spearman"
            else:
                pref = base_clean if base_clean else "var"

            if suf_clean:
                nuevo = f"{pref}_{suf_clean}"
            else:
                nuevo = pref

            nuevos_nombres.append(nuevo)

        df.columns = nuevos_nombres
        # Eliminamos la fila de sufijos (la primera)
        df = df.iloc[1:, :]

        # A veces hay una fila llamada "herramienta" que solo actúa de etiqueta
        if "herramienta" in df.index:
            df = df.drop(index="herramienta")

        return df

    # Si llegamos aquí, formato inesperado -> mensaje de ayuda
    raise SystemExit(
        "ERROR: No se pudo interpretar resumen_correlacion.tsv.\n"
        f"Columnas detectadas: {df.columns.tolist()}\n"
        "Revisa el formato o envíame un fragmento actualizado del fichero."
    )


def main(carpeta_correlaciones):
    # ---- 1. Comprobar carpeta ----
    if not os.path.isdir(carpeta_correlaciones):
        raise SystemExit(f"ERROR: no existe la carpeta {carpeta_correlaciones}")

    print(f"Usando carpeta de correlaciones: {carpeta_correlaciones}")

    # Rutas de ficheros de entrada
    f_correlaciones = os.path.join(carpeta_correlaciones, "correlaciones_por_muestra.tsv")
    f_resumen = os.path.join(carpeta_correlaciones, "resumen_correlacion.tsv")
    f_wilcoxon = os.path.join(carpeta_correlaciones, "wilcoxon_correlacion.tsv")

    # ---- 2. Leer datos ----
    print("Leyendo correlaciones por muestra...")
    df_corr = pd.read_csv(f_correlaciones, sep="\t")

    print("Preparando resumen de correlaciones...")
    df_resumen = preparar_resumen(f_resumen)

    print("Leyendo Wilcoxon por pares...")
    df_wilcox = pd.read_csv(f_wilcoxon, sep="\t")

    # Nos aseguramos de que usamos sólo Pearson
    if "pearson_r" not in df_corr.columns:
        raise SystemExit("ERROR: no se encuentra la columna 'pearson_r' en correlaciones_por_muestra.tsv")

    # ---- 3. HEATMAP 1: Pearson por muestra (herramienta × muestra) ----
    print("Generando heatmap de Pearson (herramientas × muestras)...")

    # Pivotar: filas = herramienta, columnas = muestra, valores = pearson_r
    mat_pearson = df_corr.pivot(index="herramienta", columns="muestra", values="pearson_r")

    # Orden opcional: herramientas y muestras alfabéticamente
    mat_pearson = mat_pearson.sort_index().sort_index(axis=1)

    fig, ax = plt.subplots(figsize=(20, 6))

    im = ax.imshow(
        mat_pearson.values,
        aspect="auto",
        interpolation="nearest",
        vmin=0.0,
        vmax=1.0,
    )

    # Ejes
    ax.set_yticks(np.arange(mat_pearson.shape[0]))
    ax.set_yticklabels(mat_pearson.index)

    ax.set_xticks(np.arange(mat_pearson.shape[1]))
    ax.set_xticklabels(mat_pearson.columns, rotation=90, fontsize=6)

    ax.set_xlabel("Muestra")
    ax.set_ylabel("Herramienta")
    ax.set_title("Correlación Pearson por muestra (herramienta vs Gold Standard)")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pearson r (herramienta vs GS)")

    fig.tight_layout()

    out_png_pearson = os.path.join(carpeta_correlaciones, "heatmap_pearson_por_muestra.png")
    fig.savefig(out_png_pearson, dpi=300)
    plt.close(fig)
    print(f"✅ Heatmap Pearson guardado en: {out_png_pearson}")

    # ---- 4. HEATMAP 2: Wilcoxon FDR (matriz herramienta × herramienta) ----
    print("Generando heatmap de Wilcoxon (Pearson, FDR)...")

    # Nos quedamos sólo con las filas de Pearson
    df_wilcox_p = df_wilcox[df_wilcox["metrica"] == "pearson_r"].copy()

    # Sacar lista de herramientas (del resumen ya limpio)
    herramientas = df_resumen.index.tolist()

    # A estas alturas, df_resumen debe tener pearson_median
    if "pearson_median" not in df_resumen.columns:
        raise SystemExit(
            "ERROR: después de procesar resumen_correlacion.tsv sigue sin existir 'pearson_median'.\n"
            f"Columnas actuales: {df_resumen.columns.tolist()}"
        )

    median_pearson = df_resumen["pearson_median"].to_dict()

    # Matriz vacía de comparación: filas = h1, columnas = h2
    #   +1  si h1 > h2 (mediana Pearson) y p_fdr < 0.05
    #   -1  si h1 < h2 (mediana Pearson) y p_fdr < 0.05
    #    0  si no significativo
    mat_comp = pd.DataFrame(
        data=0.0,
        index=herramientas,
        columns=herramientas
    )

    alpha = 0.05
    for _, row in df_wilcox_p.iterrows():
        h1 = row["h1"]
        h2 = row["h2"]
        p_fdr = row["p_fdr"]

        if h1 not in median_pearson or h2 not in median_pearson:
            continue

        if pd.isna(p_fdr) or p_fdr >= alpha:
            continue

        m1 = median_pearson[h1]
        m2 = median_pearson[h2]

        if m1 > m2:
            mat_comp.loc[h1, h2] = 1.0
            mat_comp.loc[h2, h1] = -1.0
        elif m1 < m2:
            mat_comp.loc[h1, h2] = -1.0
            mat_comp.loc[h2, h1] = 1.0
        else:
            # empate → 0
            continue

    # Diagonal = 0
    np.fill_diagonal(mat_comp.values, 0.0)

    fig2, ax2 = plt.subplots(figsize=(7, 6))

    cmap = plt.cm.get_cmap("bwr")  # azul (-1) blanco (0) rojo (+1)

    im2 = ax2.imshow(
        mat_comp.values,
        vmin=-1,
        vmax=1,
        cmap=cmap
    )

    ax2.set_xticks(np.arange(mat_comp.shape[1]))
    ax2.set_xticklabels(mat_comp.columns, rotation=90, fontsize=7)
    ax2.set_yticks(np.arange(mat_comp.shape[0]))
    ax2.set_yticklabels(mat_comp.index, fontsize=7)

    ax2.set_xlabel("h2")
    ax2.set_ylabel("h1")
    ax2.set_title("Comparación entre herramientas (Pearson, Wilcoxon pareado, FDR<0.05)")

    # Anotar celdas con -1, 0, 1
    for i in range(mat_comp.shape[0]):
        for j in range(mat_comp.shape[1]):
            val = mat_comp.iloc[i, j]
            if val > 0:
                txt = "1"
            elif val < 0:
                txt = "-1"
            else:
                txt = "0"
            ax2.text(j, i, txt, ha="center", va="center", fontsize=8)

    cbar2 = fig2.colorbar(im2, ax=ax2)
    cbar2.set_label("h1 vs h2 (1: h1>h2, -1: h1<h2, 0: ns)")

    fig2.tight_layout()

    out_png_wilcox = os.path.join(carpeta_correlaciones, "heatmap_wilcoxon_pearson_FDR.png")
    fig2.savefig(out_png_wilcox, dpi=300)
    plt.close(fig2)
    print(f"✅ Heatmap Wilcoxon guardado en: {out_png_wilcox}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python plot_heatmaps_correlacion.py RUTA_A_correlacion_perfiles")
        print("Ejemplo:")
        print("   python plot_heatmaps_correlacion.py correlacion_perfiles")
        raise SystemExit(1)

    carpeta = sys.argv[1]
    main(carpeta)

