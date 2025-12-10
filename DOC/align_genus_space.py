#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import sys
import pandas as pd

FIXED_COLS = {"genus", "sum_abundance"}

def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Alinea varias tablas de abundancias por género al MISMO espacio de géneros.\n"
            "No elimina géneros: usa la UNIÓN de géneros de todos los archivos.\n"
            "Cada archivo debe tener columnas: 'genus', 'sum_abundance' y columnas de muestras."
        )
    )
    ap.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Lista de TSV de entrada (GS + herramientas)."
    )
    ap.add_argument(
        "--outdir",
        required=True,
        help="Directorio donde escribir las tablas alineadas."
    )
    ap.add_argument(
        "--suffix",
        default=".aligned.tsv",
        help="Sufijo para las salidas (por defecto: .aligned.tsv)."
    )
    return ap.parse_args()

def read_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t", dtype=float)
    except Exception:
        # si falla como float, reintenta como str y luego convierte
        df = pd.read_csv(path, sep="\t", dtype=str)
        # convertir numéricos después
        for c in df.columns:
            if c not in FIXED_COLS:
                df[c] = pd.to_numeric(df[c], errors="coerce")

    cols = set(df.columns)
    if "genus" not in cols:
        sys.exit(f"[ERROR] {path}: falta columna 'genus'")
    if "sum_abundance" not in cols:
        print(f"[WARN] {path}: no tiene 'sum_abundance'; la recalcularé.", file=sys.stderr)

    # columnas de muestra = todo menos genus y sum_abundance
    sample_cols = [c for c in df.columns if c not in FIXED_COLS]

    if not sample_cols:
        sys.exit(f"[ERROR] {path}: no hay columnas de muestra (solo genus/sum_abundance).")

    # Asegurar tipo numérico en muestras
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    # Asegurar que genus es string “limpio”
    df["genus"] = df["genus"].astype(str).str.strip()

    return df, sample_cols

def main():
    args = parse_args()

    in_paths = [Path(p) for p in args.inputs]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Leer todas las tablas y recolectar géneros + columnas de muestra
    dfs = {}
    sample_cols_dict = {}
    all_genera = set()

    for p in in_paths:
        df, sample_cols = read_table(p)
        dfs[p] = df
        sample_cols_dict[p] = sample_cols
        all_genera.update(df["genus"].unique())
        print(f"[INFO] {p.name}: {len(df)} géneros, {len(sample_cols)} muestras.", file=sys.stderr)

    all_genera = sorted(all_genera)
    print(f"[INFO] Géneros totales (unión de todos los archivos): {len(all_genera)}", file=sys.stderr)

    # 2) Reindexar cada tabla al espacio global de géneros
    for p in in_paths:
        df = dfs[p]
        sample_cols = sample_cols_dict[p]

        # usar genus como índice
        df = df.set_index("genus")

        # reindex al conjunto global de géneros; géneros ausentes -> NaN
        df = df.reindex(all_genera)

        # Rellenar NaN en muestras con 0
        for c in sample_cols:
            df[c] = df[c].fillna(0.0)

        # Recalcular sum_abundance como suma de muestras
        df["sum_abundance"] = df[sample_cols].sum(axis=1)

        # Volver a tener 'genus' como columna
        df = df.reset_index()
        df.rename(columns={"index": "genus"}, inplace=True)

        # Ordenar columnas: genus, sum_abundance, luego muestras (orden alfabético)
        sample_cols_sorted = sorted(sample_cols)
        df = df[["genus", "sum_abundance"] + sample_cols_sorted]

        # Escribir salida
        out_name = p.stem + args.suffix
        out_path = outdir / out_name
        df.to_csv(out_path, sep="\t", index=False, float_format="%.16g")
        print(f"[OK] {p.name} -> {out_path.name} ({len(df)} géneros)", file=sys.stderr)


if __name__ == "__main__":
    main()