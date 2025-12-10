#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import argparse
import pandas as pd
from pathlib import Path

COL_GENUS = "genus"
COL_SUM = "sum_abundance"


def guess_is_header_line(line: str) -> bool:
    parts = line.rstrip("\n").split("\t")
    if not parts:
        return False
    if parts[0] == "" and any(".tre" in p for p in parts[1:]):
        return True
    if parts[0].lower() == "genus":
        return True
    return False


def clean_sample(filename: str) -> str:
    """
    De 'ALB.17_ERR1713331_db11.ganon_report.tre' -> 'ALB.17'
    De 'DNK.71_RA_ERR1713349_db11.ganon_report.tre' -> 'DNK.71_RA'
    """
    base = Path(filename).name
    base = re.sub(r'\.ganon_report\.tre$', '', base)
    base = re.sub(r'_ERR[^_]*.*$', '', base)
    base = re.sub(r'_db\d+$', '', base)
    return base


# -------------------------------
#   TAXDUMP: FILTRAR SOLO BACTERIAS
# -------------------------------

def load_nodes(nodes_path: Path):
    """
    Carga nodes.dmp -> diccionarios:
      - parents[taxid] = parent_taxid
      - ranks[taxid] = rank (e.g. 'genus', 'species'...)
    """
    parents = {}
    ranks = {}
    with nodes_path.open() as f:
        for line in f:
            # Formato: tax_id | parent_tax_id | rank | ...
            parts = line.split("|")
            taxid = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank = parts[2].strip()
            parents[taxid] = parent
            ranks[taxid] = rank
    return parents, ranks


def build_is_bacteria_checker(parents: dict):
    """
    Devuelve una función is_bacteria(taxid) que comprueba si un taxid
    desciende de Bacteria (taxid=2) en el árbol NCBI.
    """
    BACTERIA_TAXID = 2
    cache = {}

    def is_bacteria(tid: int) -> bool:
        if tid in cache:
            return cache[tid]
        t = tid
        # recorre hacia arriba hasta la raíz
        while True:
            if t == BACTERIA_TAXID:
                cache[tid] = True
                return True
            pt = parents.get(t)
            if pt is None or pt == t or pt in (0, 1):
                cache[tid] = False
                return False
            t = pt

    return is_bacteria


def load_bacterial_genus_names(taxdump_dir: Path, genus_names):
    """
    A partir de un directorio con taxdump (nodes.dmp, names.dmp),
    devuelve un set con los nombres (lowercase) de aquellos géneros
    que sean:
      - rank == 'genus'
      - dentro de Bacteria
      - y cuyo nombre científico coincide con los nombres de 'genus_names'
    """
    nodes_path = taxdump_dir / "nodes.dmp"
    names_path = taxdump_dir / "names.dmp"

    if not nodes_path.exists() or not names_path.exists():
        sys.exit(f"No encuentro 'nodes.dmp' o 'names.dmp' en {taxdump_dir}")

    parents, ranks = load_nodes(nodes_path)
    is_bacteria = build_is_bacteria_checker(parents)

    # Trabajamos en minúsculas para comparar nombres
    genus_lower = {g.lower() for g in genus_names}
    valid_bact_genus = set()

    with names_path.open() as f:
        for line in f:
            # names.dmp:
            # tax_id | name_txt | unique name | name class |
            parts = line.split("|")
            taxid = int(parts[0].strip())
            name_txt = parts[1].strip()
            name_class = parts[3].strip()

            if name_class != "scientific name":
                continue

            nlow = name_txt.lower()
            if nlow not in genus_lower:
                continue

            if ranks.get(taxid) != "genus":
                continue

            if not is_bacteria(taxid):
                continue

            valid_bact_genus.add(nlow)

    return valid_bact_genus


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Convierte tabla de conteos por género (ganon) a abundancias relativas por muestra + sum_abundance.\n"
            "Opcionalmente, filtra para quedarse solo con géneros bacterianos reales usando un taxdump de NCBI."
        )
    )
    ap.add_argument("input", help="Fichero TSV de entrada (usa '-' para stdin).")
    ap.add_argument(
        "-o", "--output",
        default="abundance.tsv",
        help="TSV de salida (por defecto: abundance.tsv)"
    )
    ap.add_argument(
        "--drop-root",
        action="store_true",
        help="Eliminar la fila 'root' si existe."
    )
    ap.add_argument(
        "--taxdump-dir",
        metavar="DIR",
        help=(
            "Directorio con taxdump de NCBI (names.dmp, nodes.dmp, ...). "
            "Si se proporciona, se filtra para quedarse solo con géneros bacterianos 'reales'."
        )
    )
    args = ap.parse_args()

    # Leer entrada
    raw = sys.stdin.read().strip("\n") if args.input == "-" else Path(args.input).read_text()
    lines = [ln for ln in raw.splitlines() if ln.strip() != ""]
    if not lines:
        sys.exit("Entrada vacía.")

    # Cabecera
    if guess_is_header_line(lines[0]):
        header = lines[0].split("\t")
        data_lines = lines[1:]
    else:
        first_data = lines[0].split("\t")
        n_cols = len(first_data)
        header = [COL_GENUS] + [f"sample_{i}" for i in range(1, n_cols)]
        data_lines = lines

    # Limpiar nombres de muestra
    if header[0] == "":
        header[0] = COL_GENUS
    samples = header[1:]
    clean_samples = [clean_sample(s) for s in samples]

    # DataFrame
    rows = []
    for ln in data_lines:
        parts = ln.rstrip("\n").split("\t")
        if len(parts) != len(header):
            if len(parts) < len(header):
                parts = parts + ["0"] * (len(header) - len(parts))
            else:
                parts = parts[:len(header)]
        rows.append(parts)

    df = pd.DataFrame(rows, columns=[COL_GENUS] + clean_samples)

    # Numérico
    sample_cols = [c for c in df.columns if c != COL_GENUS]
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    # Quitar 'root' si procede
    if args.drop_root:
        df = df[df[COL_GENUS] != "root"].copy()

    # --- FILTRAR SOLO GÉNEROS BACTERIANOS (SI HAY TAXDUMP) ---
    if args.taxdump_dir:
        taxdump_path = Path(args.taxdump_dir)
        valid_genus_lower = load_bacterial_genus_names(
            taxdump_path,
            df[COL_GENUS].unique()
        )
        df = df[df[COL_GENUS].str.lower().isin(valid_genus_lower)].copy()

    # --- EXCLUIR Candidatus ANTES DE ABUNDANCIAS ---
    # Eliminamos filas cuyo género empiece por "Candidatus" (mayúsc/minus da igual)
    genus_series = df[COL_GENUS].astype(str).str.strip()
    mask_candidatus = genus_series.str.contains(r"(?i)^candidatus\b", na=False)
    if mask_candidatus.any():
        n_cand = mask_candidatus.sum()
        print(
            f"[INFO] Excluyo {n_cand} filas con género que empieza por 'Candidatus'.",
            file=sys.stderr,
        )
        df = df[~mask_candidatus].copy()

    # Abundancias relativas por muestra
    col_sums = df[sample_cols].sum(axis=0)
    rel = df[sample_cols].div(col_sums.replace({0: 1}), axis=1)

    # sum_abundance = SUMA de abundancias relativas del género en TODAS las muestras
    df[COL_SUM] = rel.sum(axis=1)

    # Salida final
    out = pd.concat([df[[COL_GENUS, COL_SUM]], rel], axis=1)
    out = out.sort_values(COL_GENUS, ascending=True).reset_index(drop=True)

    # Ordenar columnas: primero 'genus', luego 'sum_abundance'
    # y después las columnas de muestras en orden alfabético
    fixed_cols = ["genus", "sum_abundance"]
    sample_cols_sorted = sorted([c for c in out.columns if c not in fixed_cols])
    out = out[fixed_cols + sample_cols_sorted]

    out.to_csv(args.output, sep="\t", index=False, float_format="%.16g")


if __name__ == "__main__":
    main()

