# ==============================================================================
# UpSet plot (ComplexUpset)
# ==============================================================================

# 1) LIBRERÍAS
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("ComplexUpset", quietly = TRUE)) install.packages("ComplexUpset")

library(tidyverse)
library(ComplexUpset)

# 2) ARCHIVOS
files_list <- list(
  "Gold Standard" = "gold_standard.tsv",
  "MetaPhlAn"     = "metaphlan_db4.tsv.genus.tsv.study.tsv",
  "mOTUs"         = "motus_db8.tsv.genus.tsv.study.tsv",
  "Kraken2"       = "kraken2_db1.tsv.genus.tsv.study.tsv",
  "Bracken"       = "bracken_db2.tsv.genus.tsv.study.tsv",
  "Ganon"         = "ganon_db11.tsv.genus.tsv.study.tsv",
  "Centrifuge"    = "centrifuge_db3.tsv.genus.tsv.study.tsv",
  "Kaiju"         = "kaiju_db5.tsv.genus.tsv.study.tsv"
)

# 3) PALETA SET SIZES
colores_tfm <- c(
  "Gold Standard" = "#D3D3D3",
  "Kaiju"         = "#F0B8E0",
  "Kraken2"       = "#F0B8E0",
  "Ganon"         = "#D8D890",
  "Bracken"       = "#D8D890",
  "MetaPhlAn"     = "#88D8E0",
  "mOTUs"         = "#88B8D8",
  "Centrifuge"    = "#58D0D8"
)

# 4) COMPROBAR ARCHIVOS
paths <- unlist(files_list, use.names = FALSE)
missing_files <- paths[!file.exists(paths)]
if (length(missing_files) > 0) {
  stop(paste(
    "Faltan estos archivos:\n",
    paste(missing_files, collapse = "\n")
  ))
}

# 5) LEER GÉNEROS
lista_generos <- lapply(files_list, function(path) {
  df <- read.delim(path, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE)
  if (!"genus" %in% names(df)) {
    stop(paste("El archivo", path, "no tiene la columna 'genus'."))
  }
  unique(df$genus[!is.na(df$genus) & df$genus != ""])
})

# 6) MATRIZ BINARIA
herramientas <- names(files_list)
todos_los_generos <- unique(unlist(lista_generos))

data_upset <- data.frame(genus = todos_los_generos, stringsAsFactors = FALSE)
for (tool in herramientas) {
  data_upset[[tool]] <- data_upset$genus %in% lista_generos[[tool]]
}

# 7) TOTALES Y ORDEN
totales <- colSums(data_upset[, herramientas])
herramientas <- names(sort(totales, decreasing = TRUE))
totales <- totales[herramientas]

# Colores mapeados a group (1..n)
colores_por_group <- setNames(
  colores_tfm[herramientas],
  as.character(seq_along(herramientas))
)


# 8) UPSET PLOT

p <- upset(
  data_upset,
  herramientas,
  name = "Herramientas de Clasificación",
  width_ratio = 0.25,
  min_size = 10,
  
  
  # Barras de intersección
  
  base_annotations = list(
    "Tamaño de la intersección" = intersection_size(
      counts = TRUE,
      mapping = aes(fill = after_stat(y))
    ) +
      scale_fill_gradient(
        low  = "#E8F6F0",  # verde pastel muy claro
        high = "#6FBF9A"   # verde pastel más oscuro
      ) +
      ylab("Géneros compartidos") +
      theme(
        panel.grid.major.x = element_blank(),
        legend.position = "none"
      )
  ),
  
  
  # Set sizes
  
  set_sizes = (
    upset_set_size(
      geom = geom_bar(
        aes(fill = factor(group)),
        width = 0.65
      )
    ) +
      geom_text(
        aes(label = after_stat(count)),
        hjust = -0.15,
        size = 3.5,
        stat = "count"
      ) +
      scale_fill_manual(values = colores_por_group, guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      xlab("Total de géneros detectados") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(hjust = 1)
      )
  ),
  
  
  # Matriz
  
  matrix = intersection_matrix(
    geom = geom_point(size = 4),
    segment = geom_segment(linewidth = 1)
  ),
  
  sort_intersections = "descending",
  sort_intersections_by = "degree",
  sort_sets = FALSE
)

# 9) GUARDAR

out_file <- "UpSet_Plot_TFM_CPU_palette_with_totals.png"
png(out_file, width = 1800, height = 1100, res = 150)
print(p)
dev.off()

cat("UpSet plot generado con éxito en:", out_file, "\n")
