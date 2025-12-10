import pandas as pd
import numpy as np
import glob
import os

# Carpeta TSV
input_folder = "aligned_tables/"

# Patron para buscar todos los tsv
pattern = os.path.join(input_folder, "*.tsv")

# Recorremos todos los archivos .tsv
for filepath in glob.glob(pattern):
    filename = os.path.basename(filepath)
    print(f"\n▶ Archivo: {filename}")

    # Cargar tabla (suponiendo separador tabulador)
    df = pd.read_csv(filepath, sep="\t", index_col=0)

    # Excluir la primera columna (abundancia total acumulada)
    samples = df.iloc[:, 1:]

    # Sumar por columnas (muestras)
    col_sums = samples.sum(axis=0)

    # Comprobacion: ¿todas las muestras suman aproximadamente 1?
    all_close = np.allclose(col_sums, 1.0, atol=1e-3)

    print(f"  Todas las muestras suman ~1? {all_close}")
    print(f"  Suma minima: {col_sums.min():.4f}")
    print(f"  Suma maxima: {col_sums.max():.4f}")
