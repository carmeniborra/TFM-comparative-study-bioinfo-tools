# SCIPT PARA GENERAR EL BOXPLOT DE MÉTRICAS AGRUPADO 
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Cargar los datos
df <- read.table("metricas_por_muestra.tsv", header = TRUE, sep = "\t")

# Limpiar los nombres de las herramientas ---
df <- df %>%
  mutate(herramienta = gsub("_db.*", "", herramienta))


# 2. Normalización y transformación a formato largo
df_plot <- df %>%
  mutate(across(c(Aitchison, JSD, BrayCurtis, Spearman, F1), 
                ~ ( . - min(., na.rm = TRUE) ) / ( max(., na.rm = TRUE) - min(., na.rm = TRUE) ))) %>%
  pivot_longer(
    cols = c(Aitchison, JSD, BrayCurtis, Spearman, F1),
    names_to = "Metrica",
    values_to = "Valor_Normalizado"
  )

# 3. Definir tus colores (5 únicos para las 5 métricas)
colores_tfm <- c(
  "Aitchison"  = "#F0B8E0", 
  "JSD"        = "#D8D890", 
  "BrayCurtis" = "#88D8E0", 
  "Spearman"   = "#88B8D8", 
  "F1"         = "#58D0D8"
)

# 4. Crear el gráfico con ajustes de legibilidad para tamaño pequeño (13x8 cm)
p <- ggplot(df_plot, aes(x = herramienta, y = Valor_Normalizado, fill = Metrica)) +
  # Reducimos el tamaño de outliers y el grosor de las líneas
  geom_boxplot(outlier.size = 0.3, lwd = 0.3, alpha = 0.9, 
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colores_tfm) +
  theme_minimal(base_size = 8) + 
  labs(
    x = "Herramientas",
    y = "Rendimiento Relativo (0-1)",
    fill = "Métrica"
  ) +
  theme(
    # Rotación de etiquetas del eje X para evitar solapamiento
    axis.text.x = element_text(angle = 35, hjust = 1, size = 7), 
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(t = -5),
    panel.grid.major.x = element_line(color = "gray95"),
    plot.margin = margin(5, 5, 5, 5)
  )

# 5. Guardado con las dimensiones exactas solicitadas
ggsave("boxplot_tfm_13x8_limpio.png", plot = p, ) 
