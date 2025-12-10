#!/usr/bin/env python3
# Script para extraer la hoja de abundancias relativas del material suplementario
import pandas as pd
import sys

xlsx = sys.argv[1]          # input
tsv_out = sys.argv[2]       # output

df = pd.read_excel(xlsx, sheet_name=1)
df.to_csv(tsv_out, sep="\t", index=False)

print(f"[OK] Exportado gold standard desde la segunda hoja → {tsv_out}")