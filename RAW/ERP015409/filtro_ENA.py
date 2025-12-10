import pandas as pd
import re

tsv_path = "ENA_ERP015409.tsv"       # TSV con 'library_name' y 'run_accession'
csv_ids_path = "sra_ids_ERP015409.csv"     # CSV sin cabecera (una ID por línea)
out_tsv = "ERP015409_filtrado.tsv"
out_txt = "run_accessions.txt"

# 1) IDs
ids = (
    pd.read_csv(csv_ids_path, header=None, dtype=str)[0]
      .dropna().astype(str).str.strip()
)
ids = ids[ids.ne("")].drop_duplicates()

# 2) Patrón palabra completa (case-insensitive; quita IGNORECASE si no lo quieres)
pattern = re.compile(r"(?<!\w)(?:%s)(?!\w)" % "|".join(map(re.escape, ids)), flags=re.IGNORECASE)

# 3) Leer y filtrar
df = pd.read_csv(tsv_path, sep="\t", dtype=str)
mask = df["library_name"].astype(str).str.contains(pattern, na=False)
df_filtrado = df[mask].copy()

# 4) Guardar TSV filtrado
df_filtrado.to_csv(out_tsv, sep="\t", index=False)

# 5) Extraer run_accession (por nombre si existe; si no, primera columna)
col_ra = "run_accession" if "run_accession" in df_filtrado.columns else df_filtrado.columns[0]
(
    df_filtrado[col_ra]
      .dropna().astype(str).str.strip()
      .drop_duplicates()
      .to_csv(out_txt, index=False, header=False)
)

print(
    f"Filas originales: {len(df)} | Filtradas: {len(df_filtrado)} | "
    f"Run accessions únicos: {df_filtrado[col_ra].nunique()} | "
    f"Guardados: {out_tsv}, {out_txt}"
)
