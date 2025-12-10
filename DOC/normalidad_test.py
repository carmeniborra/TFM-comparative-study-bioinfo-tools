#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import shapiro


def main(metrics_dir: str):
    metrics_dir = Path(metrics_dir).resolve()
    metrics_file = metrics_dir / "metricas_por_muestra.tsv"

    if not metrics_file.exists():
        print(f"ERROR: No encuentro {metrics_file}", file=sys.stderr)
        sys.exit(1)

    print(f"Leyendo métricas por muestra desde: {metrics_file}")

    df = pd.read_csv(metrics_file, sep="\t")

    # Comprobaciones mínimas
    if "herramienta" not in df.columns or "muestra" not in df.columns:
        print(
            "ERROR: La tabla debe contener las columnas 'herramienta' y 'muestra'.",
            file=sys.stderr,
        )
        print(f"Columnas detectadas: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Detectar automáticamente qué columnas son métricas:
    # todas menos herramienta y muestra
    metric_cols = [c for c in df.columns if c not in ["herramienta", "muestra"]]

    if not metric_cols:
        print(
            "ERROR: No se han encontrado columnas de métricas (solo 'herramienta' y 'muestra').",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Voy a hacer Shapiro-Wilk para estas métricas: {metric_cols}")

    resultados = []

    for metrica in metric_cols:
        # Asegurar que es numérica
        try:
            df[metrica] = pd.to_numeric(df[metrica], errors="coerce")
        except Exception as e:
            print(f"AVISO: No puedo convertir la columna {metrica} a numérico: {e}", file=sys.stderr)
            continue

        for herramienta, sub in df.groupby("herramienta"):
            valores = sub[metrica].dropna().values

            n = len(valores)
            if n < 3:
                W = np.nan
                p = np.nan
            else:
                # Shapiro-Wilk
                W, p = shapiro(valores)

            resultados.append(
                {
                    "metrica": metrica,
                    "herramienta": herramienta,
                    "n_muestras": n,
                    "shapiro_W": W,
                    "shapiro_p": p,
                }
            )

    df_out = pd.DataFrame(resultados)

    out_file = metrics_dir / "normality_tests_Shapiro.tsv"
    df_out.to_csv(out_file, sep="\t", index=False)

    print(f"✅ Tests de normalidad (Shapiro-Wilk) guardados en: {out_file}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            "Uso: 02_normalidad_metricas.py <carpeta_con_metricas>\n"
            "La carpeta debe contener 'metricas_por_muestra.tsv'.",
            file=sys.stderr,
        )
        sys.exit(1)

    main(sys.argv[1])
