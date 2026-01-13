import pandas as pd
from functools import reduce
import os

# 1. Lista de tus archivos
archivos = [
    'diamond_db6_1.1.tsv', 'diamond_db6_1.2.tsv', 'diamond_db6_2.tsv',
    'diamond_db6_3.tsv', 'diamond_db6_4.tsv', 'diamond_db6_5.tsv',
    'diamond_db6_6.1.tsv', 'diamond_db6_6.2.tsv', 'diamond_db6_7.tsv'
]

# Columnas que sirven como clave (no se repiten)
columnas_clave = ['taxonomy_id', 'name', 'rank']

dfs = []

# 2. Leer cada archivo y guardarlo en una lista
print("Leyendo archivos...")
for f in archivos:
    if os.path.exists(f):
        df = pd.read_csv(f, sep='\t')
        dfs.append(df)
    else:
        print(f"Advertencia: No se encontró el archivo {f}")

# 3. Unificar todos los DataFrames usando un 'outer merge'
print("Unificando datos...")
df_final = reduce(lambda left, right: pd.merge(left, right, on=columnas_clave, how='outer'), dfs)

# 4. Opcional: Rellenar valores vacíos (NaN) con 0 
# (Útil si los números representan conteos)
df_final = df_final.fillna(0)

# 5. Guardar el resultado
nombre_salida = 'diamond_db6_unificado.tsv'
df_final.to_csv(nombre_salida, sep='\t', index=False)

print(f"¡Hecho! Archivo guardado como: {nombre_salida}")
print(f"Dimensiones finales: {df_final.shape}")
