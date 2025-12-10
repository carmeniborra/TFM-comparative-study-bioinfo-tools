#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import argparse
import pandas as pd
from pathlib import Path

COL_CLADE = "clade_name"
COL_GENUS = "genus"
COL_SUM = "sum_abundance"


def clean_sample(filename: str) -> str:
    """
    'AUS.18a_ERR1713333_db4.metaphlan' -> 'AUS.18a'
    """
    base = Path(filename).name
    base = re.sub(r'\.metaphlan$', '', base)
    base = re.sub(r'_ERR[^_]*.*$', '', base)   # corta desde _ERR... en adelante
    base = re.sub(r'_db\d+$', '', base)        # corta sufijo _db4 si queda
    return base


# ---------- CARGA TAXDUMP NCBI Y FILTRO DE GÉNEROS VÁLIDOS ----------

def load_nodes(nodes_path: Path):
    """
    Carga nodes.dmp en dos diccionarios:
      - parent[taxid] = parent_taxid
      - rank[taxid]   = rank (e.g. 'genus', 'species', ...)
    """
    parent = {}
    rank = {}
    with nodes_path.open() as f:
        for line in f:
            # tax_id \t| parent_tax_id \t| rank \t| ...
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            try:
                taxid = int(parts[0])
                parent_taxid = int(parts[1])
            except ValueError:
                continue
            r = parts[2]
            parent[taxid] = parent_taxid
            rank[taxid] = r
    return parent, rank


def build_bacteria_flag(parent: dict):
    """
    Devuelve una función is_bacteria(taxid) que comprueba si el taxid
    pertenece a la rama de Bacteria (taxid = 2) subiendo por los padres.
    """
    BACTERIA_TAXID = 2
    ROOT_TAXID = 1
    cache = {}

    def is_bacteria(tid: int) -> bool:
        if tid in cache:
            return cache[tid]
        seen = set()
        current = tid
        while True:
            if current == BACTERIA_TAXID:
                for s in seen:
                    cache[s] = True
                cache[tid] = True
                return True
            if current == ROOT_TAXID or current not in parent or current == 0:
                for s in seen:
                    cache[s] = False
                cache[tid] = False
                return False
            seen.add(current)
            current = parent[current]

    return is_bacteria


def load_valid_bacterial_genera(taxdump_dir: str):
    """
    A partir de un directorio con nodes.dmp y names.dmp (NCBI taxdump),
    construye un conjunto de nombres de géneros válidos en Bacteria.
    """
    taxdump_dir = Path(taxdump_dir)
    nodes_path = taxdump_dir / "nodes.dmp"
    names_path = taxdump_dir / "names.dmp"

    if not nodes_path.is_file():
        sys.exit(f"No se encuentra nodes.dmp en {taxdump_dir}")
    if not names_path.is_file():
        sys.exit(f"No se encuentra names.dmp en {taxdump_dir}")

    parent, rank = load_nodes(nodes_path)
    is_bacteria = build_bacteria_flag(parent)

    # 1) Taxids de rango 'genus' que pertenezcan a Bacteria
    genus_taxids = {tid for tid, r in rank.items()
                    if r == "genus" and is_bacteria(tid)}

    # 2) Cargar nombres científicos de esos taxids
    valid_genera = set()
    with names_path.open() as f:
        for line in f:
            # tax_id | name_txt | unique_name | name_class |
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            try:
                tid = int(parts[0])
            except ValueError:
                continue
            if tid not in genus_taxids:
                continue
            name_txt = parts[1]
            name_class = parts[3]
            if name_class == "scientific name":
                valid_genera.add(name_txt)

    return valid_genera


# ---------- EXTRACCIÓN DE GÉNERO DE METAPHLAN + FILTRO NCBI ----------

def extract_genus(clade: str, valid_genera: set):
    """
    Devuelve el nombre del género si la fila corresponde a:
      - Reino Bacteria (k__Bacteria)
      - Nivel de género (contiene g__ y NO contiene s__ ni t__)
      - Género presente en el conjunto de géneros válidos de NCBI
    """

    # 1) Filtrar solo Bacteria en la anotación de MetaPhlAn
    if not clade.startswith("k__Bacteria|") and clade != "k__Bacteria":
        return None

    # 2) Si hay especie/strain, descartamos (solo nivel género)
    if re.search(r'(\|s__|\|t__)', clade):
        return None

    # 3) Buscar el último taxón de nivel género (g__)
    parts = clade.split("|")
    g = None
    for p in reversed(parts):
        if p.startswith("g__"):
            g = p[3:]
            break

    if g is None:
        if clade.startswith("g__"):
            g = clade[3:]
        else:
            return None

    g = g.strip()
    if not g:
        return None

    # Si el nombre no está directamente en NCBI, intentar una versión simplificada
    if g not in valid_genera:
        simple = re.split(r"[ _]", g)[0]
        if simple in valid_genera:
            g = simple
        else:
            return None

    return g


# ---------- MAIN ----------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Convierte tabla MetaPhlAn a abundancias relativas por GÉNERO "
            "y añade sum_abundance, filtrando solo géneros bacterianos "
            "válidos según el taxdump de NCBI."
        )
    )
    ap.add_argument("input", help="Fichero TSV de MetaPhlAn (usa '-' para stdin).")
    ap.add_argument(
        "-o", "--output",
        default="metaphlan_genus_abundance.tsv",
        help="TSV de salida (por defecto: metaphlan_genus_abundance.tsv; usa '-' para stdout)"
    )
    ap.add_argument(
        "--keep-percent",
        action="store_true",
        help="Mantener porcentajes (sin renormalizar a 1). "
             "Si quieres renormalizar a 1, NO uses esta opción."
    )
    ap.add_argument(
        "--drop-zeros",
        action="store_true",
        help="Eliminar géneros con suma total 0."
    )
    ap.add_argument(
        "--taxdump-dir",
        required=True,
        help="Directorio con nodes.dmp y names.dmp de NCBI (p.ej. 2025-09-01_taxdump)."
    )
    args = ap.parse_args()

    # Cargar lista de géneros bacterianos válidos (NCBI)
    valid_genera = load_valid_bacterial_genera(args.taxdump_dir)

    # Leer tabla MetaPhlAn
    txt = sys.stdin.read() if args.input == "-" else Path(args.input).read_text()

    # Filtrar líneas vacías y comentarios de MetaPhlAn
    lines = [ln for ln in txt.splitlines() if ln.strip() != "" and not ln.startswith("#")]
    if not lines:
        sys.exit("Entrada vacía o solo cabeceras '#'.")
    header = lines[0].split("\t")
    if header[0] != COL_CLADE:
        sys.exit(f"La primera columna debe llamarse '{COL_CLADE}'. Encontrado: '{header[0]}'")

    # DataFrame
    rows = [ln.split("\t") for ln in lines[1:]]
    width = len(header)
    fixed_rows = [r + [""] * (width - len(r)) if len(r) < width else r[:width] for r in rows]
    df = pd.DataFrame(fixed_rows, columns=header)

    # Limpiar nombres de muestra
    raw_samples = header[1:]
    samples = [clean_sample(s) for s in raw_samples]
    df = df[[COL_CLADE] + raw_samples].copy()
    df.columns = [COL_CLADE] + samples

    # Convertir a numérico (MetaPhlAn da %)
    sample_cols = samples
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    # Filtrar género y extraer nombre con filtro NCBI
    df[COL_GENUS] = df[COL_CLADE].apply(lambda cl: extract_genus(cl, valid_genera))
    df = df[df[COL_GENUS].notna()].copy()

    # Agregar por género
    agg = df.groupby(COL_GENUS, as_index=False)[sample_cols].sum()

    # --- RENORMALIZAR ---
    if args.keep_percent:
        # Dejas los % tal cual (no renormalizas a 1)
        pass
    else:
        # 1) pasar de % a proporciones
        agg[sample_cols] = agg[sample_cols] / 100.0
        # 2) renormalizar cada muestra para que sume 1
        col_sums = agg[sample_cols].sum(axis=0)
        agg[sample_cols] = agg[sample_cols].div(col_sums.replace({0: 1}), axis=1)

    # sum_abundance = suma a lo largo de todas las muestras
    if len(sample_cols) == 0:
        agg[COL_SUM] = 0.0
    else:
        agg[COL_SUM] = agg[sample_cols].sum(axis=1)

    if args.drop_zeros:
        agg = agg[agg[COL_SUM] > 0].copy()

    # Salida final
    out = agg[[COL_GENUS, COL_SUM] + sample_cols].copy()
    out.rename(columns={COL_GENUS: "genus"}, inplace=True)
    out = out.sort_values("genus", ascending=True).reset_index(drop=True)

    # Ordenar columnas: primero 'genus', luego 'sum_abundance' y resto muestras
    fixed_cols = ["genus", "sum_abundance"]
    sample_cols_sorted = sorted([c for c in out.columns if c not in fixed_cols])
    out = out[fixed_cols + sample_cols_sorted]

    # Escribir
    if args.output == "-":
        out.to_csv(sys.stdout, sep="\t", index=False, float_format="%.16g")
    else:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_path, sep="\t", index=False, float_format="%.16g")


if __name__ == "__main__":
    main()


