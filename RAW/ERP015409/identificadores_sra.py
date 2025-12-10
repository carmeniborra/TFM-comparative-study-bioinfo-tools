#!/usr/bin/env python3
import pandas as pd
import os

# Leer la hoja 1 (segunda), sin cabecera
df = pd.read_excel("sup_data_3_ERP015409.xlsx", sheet_name=1, header=None)

# Tomar la fila 1 (1-based) -> índice 0
row = df.iloc[0, :]                   # fila completa
col = pd.DataFrame({"sra_id": row.tolist()})

# Guardar a xlsx y luego a csv quitando las dos primeras filas
col.to_excel("sra_ids_ERP015409.xlsx", index=False, header=False)
col.iloc[2:].to_csv("sra_ids_ERP015409.csv", index=False, header=False)
os.remove("sra_ids_ERP015409.xlsx")
