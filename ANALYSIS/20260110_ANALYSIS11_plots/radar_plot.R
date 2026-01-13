# SCRIPT PARA CREAR EL RADAR PLOT
# 1. Cargar librerías
library(fmsb)
library(dplyr)
library(tidyr)

# 2. Cargar y limpiar los datos
df <- read.table("metricas_por_muestra.tsv", header = TRUE, sep = "\t")

# Limpiamos los nombres de las herramientas eliminando "_db" y lo que sigue
df <- df %>%
  mutate(herramienta = gsub("_db.*", "", herramienta))

# Calcular medianas por herramienta ya limpia
medians <- df %>%
  group_by(herramienta) %>%
  summarise(across(c(Aitchison, JSD, BrayCurtis, Spearman, F1), median, na.rm = TRUE))

# Función para normalizar (1 es mejor, 0 es peor)
norm_best_1 <- function(x, invert = FALSE) {
  norm <- (x - min(x)) / (max(x) - min(x))
  if (invert) return(1 - norm)
  return(norm)
}

data_radar <- medians %>%
  mutate(
    Aitchison = norm_best_1(Aitchison, invert = TRUE),
    JSD = norm_best_1(JSD, invert = TRUE),
    BrayCurtis = norm_best_1(BrayCurtis, invert = TRUE),
    Spearman = norm_best_1(Spearman),
    F1 = norm_best_1(F1)
  )

# Formatear para fmsb
radar_final <- as.data.frame(data_radar[,-1])
rownames(radar_final) <- data_radar$herramienta
radar_final <- rbind(rep(1, 5), rep(0, 5), radar_final)

# 3. Definir colores (ajustados a los nombres ya limpios)
colores_herramientas <- c(
  "bracken"   = "#E090C0",
  "ganon"     = "#C070A0",
  "kaiju"     = "#D8D890",
  "kraken2"   = "#F0B8E0",
  "metaphlan" = "#58D0D8",
  "motus"     = "#88D8E0",
  "centrifuge"= "#F5AB62" 
)

col_vector <- colores_herramientas[rownames(radar_final)[-c(1,2)]]

# 4. Generar el gráfico (13x8 cm)
png("radar_limpio_tfm_13x8.png", width = 13, height = 8, units = "cm", res = 600)

par(mar = c(1, 1, 1, 6)) 

radarchart(radar_final,
           axistype = 0,    # axistype = 0 elimina los números de los ejes
           pcol = col_vector, 
           plwd = 2, 
           plty = 1,
           cglcol = "grey85", cglty = 1, cglwd = 0.5,
           vlcex = 0.7,     # Tamaño de etiquetas de métricas (esquinas)
           title = "Perfil Comparativo de Herramientas")

# Añadir leyenda a la derecha con nombres limpios
legend(x = 1.25, y = 1, 
       legend = rownames(radar_final)[-c(1,2)], 
       bty = "n", pch = 16, col = col_vector, 
       cex = 0.6, pt.cex = 1, title = "Herramienta")
