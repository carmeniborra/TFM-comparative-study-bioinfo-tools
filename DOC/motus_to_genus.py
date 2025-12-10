#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Colapsa tabla mOTUs/Taxpasta a nivel de género sumando las especies.

(Se descarta `taxonomy_id` porque no siempre existe TaxID de género en la entrada;
 para mapear TaxID de género, usar TaxonKit abajo.)

Uso:
  python motus_species_to_genus.py input.tsv --out filtered/motus_db8.tsv.genus.tsv

  # Opcional: filtrar solo géneros bacterianos válidos según taxdump NCBI:
  python motus_species_to_genus.py input.tsv --out filtered/motus_db8.tsv.genus.tsv \
      --taxdump-dir 2025-09-01_taxdump
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import re

META = {"taxonomy_id", "name", "rank"}


# ---------- TAXDUMP NCBI: CARGA Y FILTRO DE GÉNEROS BACTERIANOS ----------

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
            # Formato típico:
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
    Usa memoización para no recalcular todo el rato.
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


def load_valid_bacterial_genera(taxdump_dir: str) -> set[str]:
    """
    A partir de un directorio con nodes.dmp y names.dmp (NCBI taxdump),
    construye un conjunto de nombres de géneros válidos en Bacteria.
    Devuelve los nombres EXACTOS tal como aparecen en NCBI (Ej: 'Escherichia').
    """
    taxdump_dir = Path(taxdump_dir)
    nodes_path = taxdump_dir / "nodes.dmp"
    names_path = taxdump_dir / "names.dmp"

    if not nodes_path.is_file():
        sys.exit(f"[ERROR] No se encuentra nodes.dmp en {taxdump_dir}")
    if not names_path.is_file():
        sys.exit(f"[ERROR] No se encuentra names.dmp en {taxdump_dir}")

    parent, rank = load_nodes(nodes_path)
    is_bacteria = build_bacteria_flag(parent)

    # 1) Taxids de rango 'genus' que pertenezcan a Bacteria
    genus_taxids = {
        tid for tid, r in rank.items()
        if r == "genus" and is_bacteria(tid)
    }

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

    if not valid_genera:
        print("[WARN] No se encontraron géneros bacterianos válidos en el taxdump.", file=sys.stderr)
    else:
        print(f"[INFO] Géneros bacterianos válidos cargados desde taxdump: {len(valid_genera)}", file=sys.stderr)

    return valid_genera


# ---------- PARSEO DE ARGUMENTOS Y UTILIDADES ----------

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp", help="TSV de entrada (Taxpasta/mOTUs)")
    ap.add_argument("--out", default="motus_collapsed_genus.tsv",
                    help="TSV de salida (ancho)")
    ap.add_argument(
        "--taxdump-dir",
        default=None,
        help=(
            "Directorio con taxdump de NCBI (nodes.dmp, names.dmp, ...). "
            "Si se proporciona, solo se conservarán los géneros bacterianos válidos."
        ),
    )
    return ap.parse_args()


def extract_genus(name: str) -> str:
    if name is None:
        return ""
    s = str(name)
    # quitar prefijo s__ si aparece
    s = re.sub(r"^s__", "", s)
    # reemplazar underscores por espacios
    s = s.replace("_", " ").strip()
    # genus = primera palabra
    if not s:
        return ""
    return s.split()[0]


# ---------- MAIN ----------

def main():
    args = parse_args()

    try:
        df = pd.read_csv(args.inp, sep="\t", dtype=str)
    except Exception as e:
        print(f"[ERROR] No se pudo leer {args.inp}: {e}", file=sys.stderr)
        sys.exit(1)

    # Detectar columnas mínimas
    cols = set(df.columns)
    if not {"name", "rank"}.issubset(cols):
        print(f"[ERROR] La tabla debe tener columnas 'name' y 'rank'. Encontradas: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Quedarnos con species
    df_species = df[df["rank"].astype(str).str.strip().str.lower() == "species"].copy()
    if df_species.empty:
        print("[ERROR] No hay filas con rank=species.", file=sys.stderr)
        sys.exit(2)

    # Columnas de muestra = todas menos metadatos
    sample_cols = [c for c in df.columns if c not in META]
    if not sample_cols:
        print("[ERROR] No se detectaron columnas de muestra.", file=sys.stderr)
        sys.exit(3)

    # Extraer genus
    df_species["genus"] = df_species["name"].apply(extract_genus)
    df_species = df_species[df_species["genus"].astype(str).str.strip() != ""]
    if df_species.empty:
        print("[ERROR] No se pudieron extraer nombres de género.", file=sys.stderr)
        sys.exit(4)

    # Convertir a numérico y sumar por género
    tmp = df_species.copy()
    for c in sample_cols:
        tmp[c] = pd.to_numeric(tmp[c], errors="coerce").fillna(0)

    agg = tmp.groupby("genus", as_index=False)[sample_cols].sum()
    # ordenar por nombre de género
    agg = agg.sort_values("genus").reset_index(drop=True)

    # Construir salida ancha (name=genus, rank=genus)
    wide = agg.rename(columns={"genus": "name"})
    wide.insert(1, "rank", "genus")

    # --- FILTRAR POR TAXDUMP NCBI (OPCIONAL) ---
    if args.taxdump_dir:
        valid_genera = load_valid_bacterial_genera(args.taxdump_dir)

        g = wide["name"].astype(str).str.strip()

        # intento directo
        mask_valid = g.isin(valid_genera)

        # intento simplificado (por si vinieran sufijos raros con espacios/underscores)
        simple = g.str.split(r"[ _]", n=1).str[0]
        mask_valid |= simple.isin(valid_genera)

        before = len(wide)
        wide = wide.loc[mask_valid].copy()
        after = len(wide)
        print(
            f"[INFO] Filtrados {before - after} géneros no bacterianos o no válidos según NCBI "
            f"({after} géneros restantes).",
            file=sys.stderr,
        )

        if wide.empty:
            print("[ERROR] Tras aplicar filtro NCBI no queda ningún género válido.", file=sys.stderr)
            sys.exit(5)

    # Guardar salida ancha
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] Escrito TSV colapsado a género: {out_path} ({len(wide)} géneros)")

if __name__ == "__main__":
    main()
