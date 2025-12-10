#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch: compute per-sample relative abundances for all Taxpasta-style TSVs in a directory,
and optionally export a "Relative abundances"-style wide table per file:
columns = ['genus', 'sum_abundance', <sample1>, <sample2>, ...].

Usage (typical):
  python re_abundance_study.py /path/to/tsvs \
    --pattern "*.tsv" --outdir rel_out \
    --level genus --study-wide \
    --sample-trim "_ERR"
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import re

META = ["taxonomy_id", "name", "rank"]


def compute_relative(df, level=None, sample_trim=None):
    # Optionally filter by rank
    if level:
        df = df[df["rank"].astype(str).str.strip().str.lower() == level.strip().lower()].copy()
        if df.empty:
            return None, [], "No rows after rank filter"

    # Detect sample columns
    sample_cols = [c for c in df.columns if c not in META]
    if not sample_cols:
        return None, [], "No sample columns"

    # Optional: trim sample names
    if sample_trim:
        new_names = []
        for s in sample_cols:
            if sample_trim == "_ERR":
                # Caso especial: cortar justo antes de "_ERR"
                if "_ERR" in s:
                    left = s.split("_ERR", 1)[0]
                else:
                    left = s
            else:
                # Comportamiento genérico: cortar por el delimitador dado
                if sample_trim in s:
                    left = s.split(sample_trim, 1)[0]
                else:
                    left = s

            # Por seguridad, limpia posibles guiones/puntos/underscores al final
            left = re.sub(r"[_\.\-]+$", "", left)
            new_names.append(left)

        # Resolver duplicados añadiendo sufijos __2, __3, ...
        seen = {}
        resolved = []
        for n in new_names:
            if n in seen:
                seen[n] += 1
                resolved.append(f"{n}__{seen[n]}")
            else:
                seen[n] = 1
                resolved.append(n)

        rename_map = {old: new for old, new in zip(sample_cols, resolved)}
        df = df.rename(columns=rename_map)
        sample_cols = [rename_map[c] for c in sample_cols]

    # Numeric conversion and column-wise normalization
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    sums = df[sample_cols].sum(axis=0)
    for c in sample_cols:
        col_sum = float(sums[c])
        if col_sum == 0.0:
            df[c] = 0.0
        else:
            df[c] = df[c] / col_sum

    return df, sample_cols, None


def to_study_wide(rel_df, sample_cols, expect_level="genus"):
    """
    Build a DataFrame with columns: ['genus', 'sum_abundance'] + sample_cols
    - Uses 'name' as genus label.
    - If multiple rows share the same genus, sums their normalized abundances.
    - sum_abundance is the row-wise sum across all sample columns.
    - Sorts alphabetically by genus to mirror the paper sheet.
    - Also sorts sample columns alphabetically.
    """
    # Derive genus column from 'name'
    if "name" not in rel_df.columns:
        raise ValueError("Expected column 'name' to derive 'genus'")

    gdf = rel_df.copy()
    gdf["genus"] = gdf["name"].astype(str)

    # Optional sanity: if rank exists and expect_level is given, warn if mismatched rows exist
    if expect_level and "rank" in gdf.columns:
        # (No hard fail; just ensure we grouped correctly)
        pass

    # Collapse by genus
    collapsed = gdf.groupby("genus", as_index=False)[sample_cols].sum()

    # Ordenar alfabeticamente los nombres de muestra
    sample_cols_sorted = sorted(sample_cols)

    # Compute sum_abundance across samples
    collapsed["sum_abundance"] = collapsed[sample_cols_sorted].sum(axis=1)

    # Order columns and sort
    ordered_cols = ["genus", "sum_abundance"] + sample_cols_sorted
    collapsed = collapsed[ordered_cols].sort_values("genus", kind="mergesort").reset_index(drop=True)

    return collapsed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inpdir", type=Path, help="Directory containing TSVs")
    ap.add_argument("--pattern", default="*.tsv", help="Glob pattern (default: *.tsv)")
    ap.add_argument("--outdir", default=None, help="Output directory (default: same as inpdir)")
    ap.add_argument("--level", default=None, help="Keep only rows with rank == <level> (e.g. genus/species)")
    ap.add_argument(
        "--sample-trim",
        default=None,
        help="Delimiter to trim sample names at first occurrence (e.g. '_ERR')"
    )
    ap.add_argument("--merge-long", default=None, help="Path to write a merged long-format CSV across all files")
    ap.add_argument("--tool-label", default=None, help="Tool label for long output (same label for all files)")
    ap.add_argument("--tool-from-filename", action="store_true", help="Infer tool from filename stem")
    ap.add_argument("--long-schema", choices=["paper", "generic"], default="paper",
                    help="Header style for merged long CSV (default: paper)")
    ap.add_argument("--study-wide", action="store_true",
                    help="Also write a per-file table in the paper's 'Relative abundances' schema (genus + sum_abundance + samples).")
    args = ap.parse_args()

    files = sorted(args.inpdir.glob(args.pattern))
    if not files:
        print(f"[ERROR] No files in {args.inpdir} matching {args.pattern}", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir) if args.outdir else args.inpdir
    outdir.mkdir(parents=True, exist_ok=True)

    merged_parts = []
    written = 0
    skipped = 0

    for f in files:
        try:
            df = pd.read_csv(f, sep="\t", dtype=str)
        except Exception as e:
            print(f"[WARN] {f.name}: read error: {e}; skipping", file=sys.stderr)
            skipped += 1
            continue

        for c in ["name", "rank"]:
            if c not in df.columns:
                print(f"[WARN] {f.name}: missing required column '{c}'; skipping", file=sys.stderr)
                skipped += 1
                df = None
                break
        if df is None:
            continue

        rel_df, sample_cols, err = compute_relative(
            df,
            level=args.level,
            sample_trim=args.sample_trim
        )
        if rel_df is None:
            print(f"[INFO] {f.name}: skipped ({err})", file=sys.stderr)
            skipped += 1
            continue

        # Write normalized per-file output (.rel.tsv)
        out_name = f.with_suffix(f.suffix + ".rel.tsv") if f.suffix else Path(str(f) + ".rel.tsv")
        out_path = outdir / out_name.name
        rel_df.to_csv(out_path, sep="\t", index=False)
        written += 1

        # Write study-wide (paper-like) table (.study.tsv)
        if args.study_wide:
            try:
                study_df = to_study_wide(rel_df, sample_cols, expect_level=(args.level or "genus"))
                study_name = f.with_suffix(f.suffix + ".study.tsv") if f.suffix else Path(str(f) + ".study.tsv")
                study_path = outdir / study_name.name
                study_df.to_csv(study_path, sep="\t", index=False)
            except Exception as e:
                print(f"[WARN] {f.name}: could not write study-wide table: {e}", file=sys.stderr)

        # Add to merged long output if requested
        if args.merge_long:
            long = rel_df.melt(
                id_vars=[col for col in rel_df.columns if col not in sample_cols],
                value_vars=sample_cols,
                var_name="Sample",
                value_name="Relative abundance"
            )
            # Determine tool label
            if args.tool_label:
                tool = args.tool_label
            elif args.tool_from_filename:
                tool = f.stem
            else:
                tool = "unknown"

            if args.long_schema == "paper":
                # Paper schema: Source file, Tool, Sample, Taxon, Rank, Relative abundance
                long = long.rename(columns={"name": "Taxon", "rank": "Rank"})
                long.insert(0, "Tool", tool)
                long.insert(0, "Source file", f.name)
                # Ensure ordering of columns
                cols = ["Source file", "Tool", "Sample", "Taxon", "Rank", "Relative abundance"]
                # Retain only expected columns (if taxonomy_id exists, drop it)
                long = long[[c for c in cols if c in long.columns] +
                            [c for c in long.columns if c not in cols]]
                if all(c in long.columns for c in cols):
                    long = long[cols]
            else:
                # Generic schema: source_file, tool, sample, name, rank, rel_abundance
                long.insert(0, "tool", tool)
                long.insert(0, "source_file", f.name)
                long = long.rename(columns={"Sample": "sample", "Relative abundance": "rel_abundance"})

            merged_parts.append(long)

    if args.merge_long and merged_parts:
        merged = pd.concat(merged_parts, ignore_index=True)
        Path(args.merge_long).parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(args.merge_long, index=False)
        print(f"[OK] Wrote merged long CSV: {args.merge_long} ({len(merged)} rows).")

    print(f"[DONE] Files processed: {len(files)} | written: {written} | skipped: {skipped} | outdir: {outdir}")


if __name__ == "__main__":
    main()
