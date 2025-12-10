#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filtra TSV para quedarse solo con filas donde rank == 'genus' y,
opcionalmente, solo con géneros bacterianos válidos según el taxdump de NCBI.

Uso:
  # Solo rank=genus
  python to_genus.py /ruta/a/tsvs --pattern "*.tsv" --outdir filtered

  # Rank=genus + filtro NCBI (solo géneros bacterianos válidos)
  python to_genus.py /ruta/a/tsvs --pattern "*.tsv" --outdir filtered \
      --taxdump-dir 2025-09-01_taxdump
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


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


# ---------- FILTRO POR RANK=GENUS (+ OPCIONAL NCBI) ----------

def filter_one(path: Path, valid_genera: set[str] | None = None) -> pd.DataFrame | None:
    """Lee un TSV/TSV.GZ, filtra genus y (opcionalmente) por géneros válidos de NCBI,
    y excluye géneros cuyo nombre empiece por 'Candidatus' (columna 'name')."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, compression="infer")
    except Exception as e:
        print(f"[WARN] No se pudo leer {path.name}: {e}", file=sys.stderr)
        return None

    if "rank" not in df.columns:
        print(f"[WARN] {path.name}: no existe la columna 'rank'; lo salto.", file=sys.stderr)
        return None

    # 1) rank == genus
    mask = df["rank"].astype(str).str.strip().str.lower() == "genus"
    out = df.loc[mask].copy()
    if out.empty:
        print(f"[INFO] {path.name}: 0 filas con rank=genus.", file=sys.stderr)
        return None

    # Nos aseguramos de que existe la columna 'name' (donde están los géneros)
    if "name" not in out.columns:
        print(
            f"[WARN] {path.name}: no existe columna 'name'; "
            f"no puedo aplicar filtro NCBI ni excluir Candidatus, dejo solo rank=genus.",
            file=sys.stderr,
        )
        return out

    # 2) Si se ha proporcionado taxdump, filtrar por géneros bacterianos válidos
    if valid_genera is not None:
        g = out["name"].astype(str).str.strip()

        # Intento directo
        mask_valid = g.isin(valid_genera)

        # Intento simplificado por si hay sufijos tipo 'Escherichia_unclassified'
        simple = g.str.split(r"[ _]", n=1).str[0]
        mask_valid |= simple.isin(valid_genera)

        out = out.loc[mask_valid].copy()

        if out.empty:
            print(
                f"[INFO] {path.name}: 0 filas con rank=genus y género bacteriano válido según NCBI.",
                file=sys.stderr,
            )
            return None

    # 3) EXCLUIR Candidatus usando la columna 'name'
    g2 = out["name"].astype(str).str.strip()

    # Filtra cualquier nombre que empiece por "Candidatus" (may/max da igual)
    mask_candidatus = g2.str.contains(r"(?i)^candidatus\b", na=False)

    if mask_candidatus.any():
        n_cand = mask_candidatus.sum()
        print(
            f"[INFO] {path.name}: excluyo {n_cand} filas con nombre que empieza por 'Candidatus'.",
            file=sys.stderr,
        )

    out = out.loc[~mask_candidatus].copy()

    if out.empty:
        print(
            f"[INFO] {path.name}: 0 filas tras aplicar filtros (genus + NCBI + excluir Candidatus).",
            file=sys.stderr,
        )
        return None

    return out



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inpdir", type=Path, help="Directorio con TSVs")
    ap.add_argument(
        "--pattern",
        default="*.tsv",
        help="Patrón de archivos (glob), p.ej. '*.tsv' o '*.tsv.gz'"
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help="Directorio de salida; por defecto, el mismo que inpdir"
    )
    ap.add_argument(
        "--taxdump-dir",
        default=None,
        help=(
            "Directorio con taxdump de NCBI (nodes.dmp, names.dmp, ...). "
            "Si se proporciona, solo se conservarán los géneros válidos de Bacteria."
        ),
    )
    args = ap.parse_args()

    files = sorted(args.inpdir.glob(args.pattern))
    if not files:
        print(f"[ERROR] No se encontraron archivos en {args.inpdir} con patrón {args.pattern}", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir) if args.outdir else args.inpdir
    outdir.mkdir(parents=True, exist_ok=True)

    # Cargar, si procede, lista de géneros bacterianos válidos (NCBI)
    valid_genera = None
    if args.taxdump_dir:
        valid_genera = load_valid_bacterial_genera(args.taxdump_dir)

    count_written = 0

    for f in files:
        out_df = filter_one(f, valid_genera=valid_genera)
        if out_df is None:
            continue

        # Guardar por archivo
        if f.suffix:  # .tsv o .tsv.gz (se conserva el sufijo principal)
            out_name = f.with_suffix(f.suffix + ".genus.tsv")
        else:
            out_name = Path(str(f) + ".genus.tsv")
        out_path = outdir / out_name.name
        out_df.to_csv(out_path, sep="\t", index=False)
        count_written += 1

    print(f"[OK] Archivos filtrados y guardados en {outdir}: {count_written}", file=sys.stderr)


if __name__ == "__main__":
    main()

