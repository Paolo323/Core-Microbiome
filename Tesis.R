library(dada2)
library(DECIPHER)
library(DESeq2)
library(tidyverse)
library(phyloseq)
library(ggtree)
library(ggplot2)
library(plyr)
library(metagMisc)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(vegan)
library(GGally)
library(viridis)
library(ranacapa)
library(patchwork)
library(microbiome)
library(microbenchmark)
library(ape)
library(khroma)
library(microViz)
library(ComplexHeatmap)
library(ggtext)
library(ggraph)
library(DT)
library(corncob)
library(cluster)
library(reshape2)
library(pheatmap)
library(gridExtra)
library(grid)
library(scales)
library(colorspace)
#################
# Crear la lista de objetos phyloseq
phylo_list <- list(
  Estudio_1 = psEstudio1,
  Estudio_2 = psEstudio3,
  Estudio_3 = psEstudio5,
  Estudio_4 = psEstudio6,
  Estudio_5 = psEstudio2,
  Estudio_6 = psEstudio8,
  Estudio_7 = psEstudio9
)
# Figura 1: Curvas de acumulación de ASV con profundidad de secuenciación
rarefaction_fig <- list() # Crear una lista vacía para las figuras
uniq <- names(phylo_list) # Hacer una lista de nombres de estudios

# Paleta de colores
palette <- brewer.pal(10, "Paired")
palette <- palette[c(-1, -8)] # Ajustar la paleta de colores si es necesario

# Calcular el máximo número de ASVs detectados entre todos los estudios
max_y <- max(sapply(phylo_list, function(data) max(taxa_sums(data))))

for (i in 1:length(uniq)) { # Para cada estudio i
  data_1 <- phylo_list[[i]] # Subconjunto del objeto phyloseq para el estudio i
  
  # Comprobar si los datos están vacíos o si todos los taxones tienen conteos cero
  if (is.null(data_1) || all(taxa_sums(data_1) == 0)) {
    warning(paste("El objeto phyloseq para", uniq[i], "es NULL o tiene suma de tasas igual a 0."))
    rarefaction_fig[[i]] <- NULL
    next
  }
  
  # Filtrar datos con conteos insuficientes
  min_count <- 3
  otu_table_filtered <- otu_table(data_1)[rowSums(otu_table(data_1)) >= min_count, ] 
  
  # Verificar dimensiones de la tabla filtrada
  if (nrow(otu_table_filtered) == 0) {
    warning(paste("Después de filtrar, el objeto phyloseq para", uniq[i], "tiene 0 filas."))
    rarefaction_fig[[i]] <- NULL
    next
  }
  
  # Crear un objeto phyloseq filtrado
  data_1_filtered <- phyloseq(otu_table(otu_table_filtered), sample_data(data_1), tax_table(data_1))
  
  # Crear el gráfico de rarefacción con ranacapa::ggrare
  tryCatch({
    p <- ranacapa::ggrare(data_1_filtered, step = 100, se = FALSE, plot = FALSE) +
      xlim(c(0, 15000)) # Limitar el eje x a 15,000 lecturas por muestra
    p1 <- p + geom_line(color = palette[i], size = 0.1) +
      geom_vline(xintercept = 10000) +
      geom_hline(yintercept = 200, linetype = "dashed", size = 0.5) +
      theme(legend.position = "none") +
      theme_bw(base_size = 11) +
      theme(axis.title.y = element_blank()) +
      xlab(uniq[i]) +
      theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
    rarefaction_fig[[i]] <- p1 # Guardar la figura en la lista
  }, error = function(e) {
    warning(paste("Error en el gráfico de rarefacción para", uniq[i], ":", e$message))
    rarefaction_fig[[i]] <- NULL
  })
}

# Filtrar las figuras que no son NULL y que son objetos gráficos válidos
rarefaction_fig <- rarefaction_fig[!sapply(rarefaction_fig, is.null) & 
                                     sapply(rarefaction_fig, inherits, what = "gg")]

# Verificar el número de gráficos restantes después de filtrar
num_figs <- length(rarefaction_fig)

if (num_figs > 0) {
  # Calcular el número de filas y columnas en función del número de gráficos restantes
  ncol <- 3
  nrow <- ceiling(num_figs / ncol)
  
  # Crear la figura compuesta usando ggarrange solo con las figuras válidas
  A <- ggarrange(plotlist = rarefaction_fig, 
                 ncol = ncol, nrow = nrow)
  
  # Añadir anotaciones a la figura compuesta
  A <- annotate_figure(A, bottom = text_grob("Profundidad de secuenciación", size = 14),
                       left = text_grob("ASVs detectados", rot = 90, size = 14),
                       top = text_grob("", size = 14))
  
  # Mostrar la figura
  print(A)
} else {
  message("No hay gráficos válidos para mostrar.")
}
################################################################################
###############
# Figura 2: Curvas de acumulación de ASV con tamaño de muestra
# Crear listas vacías para los resultados
SAClist <- list()
poolaccum_list <- list()
# Hacer una lista de nombres de estudios
uniq <- names(phylo_list)
for (i in 1:length(uniq)) { # Para cada estudio i 
  data_1 <- phylo_list[[i]] # Subconjunto del objeto phyloseq para el estudio i
  # Eliminar cualquier rastro de taxones que ya no estén presentes en el conjunto de datos
  data_1 <- prune_taxa(taxa_sums(data_1) > 0, data_1)
  # Transponer la tabla OTU
  data_1_matrix <- data.frame(t(otu_table(data_1)))
  # Aplicar specaccum() para calcular la acumulación de ASVs
  data_1_specaccum <- vegan::specaccum(data_1_matrix, method = "exact", permutations = 999)
  # camino principal 
  sac_df<- data_1_specaccum$sites ##sites = samples
  sac_df<-data.frame(sac_df)
  names(sac_df)[1]<-"Site"
  sac_df$Richness <- data_1_specaccum$richness #import ASV richness to dataframe
  sac_df$SD <- data_1_specaccum$sd 
  ####### Convertir el resultado en un marco de datos
  sac_df <- data.frame(Site = data_1_specaccum$sites, 
                       Richness = data_1_specaccum$richness, 
                       SD = data_1_specaccum$sd)
  # Estimar el número TOTAL de ASVs
  sac_total_estimated <- vegan::specpool(data_1_matrix)
  sac_df$Total <- sac_total_estimated$boot
  sac_df$Total_jackknife <- sac_total_estimated$jack1
  sac_df$Total_chao <- sac_total_estimated$chao
  sac_df$Study <- as.character(uniq[[i]]) # Cambiar nombre de la variable a Study
  # Añadir el marco de datos a la lista
  SAClist[[i]] <- sac_df
  # Calcular el efecto del tamaño de muestra en la estimación total
  poolaccum_df <- poolaccum(data_1_matrix, permutations = 1000, minsize = 3)
  poolaccum_df_means <- data.frame(poolaccum_df$means)
  poolaccum_df_means$Study <- uniq[[i]] # Cambiar nombre de la variable a Study
  # Añadir el marco de datos a la lista
  poolaccum_list[[i]] <- poolaccum_df_means
}
# Combinar todos los marcos de datos en uno solo
names(SAClist) <- uniq
sac_df_all <- do.call(rbind, SAClist)
poolaccum_df_all <- do.call(rbind, poolaccum_list)
# Ajustar el formato para las figuras
poolaccum_df_all$Study <- factor(poolaccum_df_all$Study, levels = uniq)
# Crear un marco de datos con los totales estimados por estudio
study_totals <- sac_df_all %>%
  distinct(Total_jackknife, .keep_all = TRUE)
study_totals$Site <- max(sac_df_all$Site) + 1 # Asignar un tamaño de muestra mayor para el total
study_totals$Richness <- NA
study_totals$SD <- NA
# Combinar los marcos de datos
sac_df_fig <- rbind(sac_df_all, study_totals)
sac_df_fig$Study <- factor(sac_df_fig$Study, levels = uniq)
# Establecer paleta de colores
palette <- brewer.pal(10, "Paired")
palette <- palette[c(-1, -8)]
# Crear la Figura 2
B <- ggplot(sac_df_fig, aes(x = Site, y = Richness, group = Study)) +
  geom_line(alpha = 1, linetype = "dashed", size = 1, aes(col = Study)) +
  geom_point(aes(col = Study), size = 1.5, alpha = 0.7) +
  theme_bw() +
  ylim(0, 500) +
  xlab("N") +
  ylab("Número de ASVs") +
  theme(text = element_text(size = 14)) +
  scale_color_manual(values = palette) +
  ggtitle("") +
  guides(colour = guide_legend(title = "Estudio", override.aes = list(size = 7)))
B <- B + scale_x_log10() + scale_y_log10() + ggtitle("A") + 
  theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 16, face = "bold"),
        legend.position = "none") # Eliminar leyenda de este gráfico
################################################################################
#############
####### Figura 3: Proporción acumulada de ASV detectado con el tamaño de la muestra
summary<-sac_df_all %>% group_by(Study) %>%
  summarize(max_asvs = max(Richness), predicted = max(Total_jackknife))
summary$percent<-summary$max_asvs/summary$predicted
sac_df_all$Percent<-sac_df_all$Richness/sac_df_all$Total_jackknife
sac_df_all$Study<-factor(sac_df_all$Study, levels = uniq)
segments_df<-subset(sac_df_all, Percent <0.51 & Percent > 0.49)
segments_df<-segments_df%>%distinct(Study, .keep_all =TRUE)
C <- ggplot(sac_df_all, aes(x = Site, y = Percent, group = Study)) +
  geom_line(alpha = 0.5, linetype = "dashed") +
  geom_point(aes(col = Study), size = 1.5, alpha = 0.9) +
  theme_bw() +
  xlab("N") +
  ylab("Proporción de ASV detectados") +
  ylim(0, 1) +  # Limitar el eje Y para mayor contraste en proporciones más bajas
  theme(text = element_text(size = 14)) +
  scale_color_manual(values = palette) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  guides(colour = guide_legend(title = "Estudio", override.aes = list(size = 7))) +
  ggtitle("") +
  geom_segment(data = segments_df, aes(x = Site, xend = Site, y = 0, yend = 0.5, col = Study), 
               size = 1.5, linetype = "dashed") +
  scale_x_log10() + ggtitle("B") + 
  theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 16, face = "bold"),
        legend.position = "bottom") # Eliminar leyenda de este gráfico
################################################################################
# Figura N°4: Efecto del tamaño de la muestra en el grupo regional de ASVs
max_estimates<-poolaccum_df_all %>% 
  group_by(Study) %>%
  summarise(max_est = max(Jackknife.1), 
            max_n = max(N))
D <- ggplot(poolaccum_df_all, 
            aes(x = N, 
                y = Jackknife.1, 
                group = Study)) +
  geom_point(aes(col = Study),
             size = 2, alpha = 0.7) +  # Aumentar tamaño de puntos
  geom_line(aes(col = Study),
            alpha = 0.5) +  # Conectar puntos con líneas
  scale_color_manual(values = palette) +
  theme_bw(base_size = 14) +
  guides(colour = guide_legend(title = "Estudio", override.aes = list(size = 7))) +
  ggtitle("") +
  ylab("N° ASVs estimados") +
  scale_x_log10() +  # Escala logarítmica en el eje X
  geom_hline(yintercept = 300, linetype = "dashed", color = "grey20") + ggtitle("C") + 
  theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 16, face = "bold"), legend.position = "none")
################################################################################
# Figura 5: Estimacion de la prevalencia y abundancia de ASV
# Función para calcular la prevalencia
prevalence <- function(physeq, add_tax = TRUE){
  trows <- taxa_are_rows(physeq)
  otutab <- as.data.frame(otu_table(physeq))
  
  # Transponer tabla OTU si los taxones no son filas
  if(trows == FALSE){
    otutab <- t(otutab)
  }
  
  # Estimar prevalencia (número de muestras con OTU presente)
  prevdf <- apply(X = otutab, MARGIN = 1, FUN = function(x){sum(x > 0)})
  
  # Agregar recuentos de abundancia total y media por OTU
  prevdf <- data.frame(Prevalence = prevdf,
                       TotalAbundance = taxa_sums(physeq),
                       MeanAbundance = rowMeans(otutab),
                       MedianAbundance = apply(otutab, 1, median))
  
  # Agregar tabla de taxonomía si está presente
  if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(physeq))
  }
  return(prevdf)
}

# Crear una lista de objetos phyloseq por estudio y calcular la prevalencia
Prevlist <- list()
uniq <- names(phylo_list)
for (i in 1:length(uniq)){
  data_1 <- phylo_list[[i]]
  data_1 <- rarefy_even_depth(data_1, sample.size = 15000, rngsee = 100, replace = TRUE, trimOTUs = FALSE, verbose = FALSE)
  occupancy_abundance <- prevalence(data_1)
  occupancy_abundance$host_study <- as.character(uniq[[i]])
  
  # Calcular abundancia relativa y prevalencia relativa
  occupancy_abundance$RelAbundance <- (occupancy_abundance$TotalAbundance / sum(occupancy_abundance$TotalAbundance))
  occupancy_abundance$RelPrev <- (occupancy_abundance$Prevalence / length(sample_data(data_1)$Run))
  occupancy_abundance$ASV <- taxa_names(data_1)
  occupancy_abundance$Sample_size <- length(unique(sample_data(data_1)$Bases))
  Prevlist[[i]] <- occupancy_abundance
}

# Unir los resultados en un solo dataframe
occupancy_abundance_df <- do.call(rbind, Prevlist)
occupancy_abundance_df$host_study <- factor(occupancy_abundance_df$host_study, levels = uniq)

# Definir umbrales de prevalencia en términos de porcentaje
prev_thresholds <- c(0.3, 0.5, 0.7, 0.9)  # 30%, 50%, 70%, 90%

# Clasificar ASVs según los umbrales de prevalencia
occupancy_abundance_df$Prevalence_alt <- cut(occupancy_abundance_df$RelPrev, 
                                             breaks = c(0, 0.3, 0.5, 0.7, 0.9, 1),
                                             labels = c("Abajo_30%", "30-50%", "50-70%", "70-90%", "Arriba_90%"))

# Calcular la proporción de ASVs en cada categoría de prevalencia por estudio
proportion_singletons <- data.frame(table(occupancy_abundance_df$Prevalence_alt, occupancy_abundance_df$host_study))
names(proportion_singletons) <- c("Prevalence", "Study", "Count")

# Remover las prevalencias con valor cero (si existen)
proportion_singletons <- subset(proportion_singletons, Prevalence != "0")

# Crear la visualización ajustada con los umbrales de prevalencia en porcentaje
E <- ggplot(proportion_singletons, aes(fill = Prevalence, y = Count, x = Study)) +
  geom_bar(position = "fill", stat = "identity") +
  ylab("Frecuencia porcentual de ASVs") +
  theme_grey(base_size = 12) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  guides(fill = guide_legend(title = "Prevalencia")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme() + ggtitle("A") + 
  theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 16, face = "bold"), legend.position = "none")
################################################
##############
###############
# Figura 6: Proporción de ASV por muestra detectada 
# Crear la lista de figuras individuales y el dataframe de prevalencia por muestra
fig_list <- list()
prevalence_barplot_df <- list()

for (i in 1:length(uniq)) {
  data_1 <- phylo_list[[i]]  # Subconjunto del objeto phyloseq para cada estudio
  data_1 <- prune_taxa(taxa_sums(data_1) > 0, data_1)
  
  # Agregar categoría de prevalencia a la tabla taxonómica como columna
  taxtable <- data.frame(tax_table(data_1))
  taxtable$ASV <- row.names(taxtable)
  
  oc_df <- subset(occupancy_abundance_df, host_study == uniq[i])
  oc_df <- oc_df[, c("ASV", "Prevalence_alt")]
  
  new_taxtable <- merge(taxtable, oc_df, by = "ASV", all.x = TRUE)
  row.names(new_taxtable) <- new_taxtable$ASV
  new_taxtable <- new_taxtable[match(row.names(taxtable), row.names(new_taxtable)), ]
  new_taxtable <- new_taxtable[, c("Kingdom", "Prevalence_alt")]
  new_taxtable <- as.matrix(new_taxtable)
  tax_table(data_1) <- tax_table(new_taxtable)
  
  # Presencia-ausencia para que la abundancia sea igual a la prevalencia
  data_pa <- phyloseq_standardize_otu_abundance(data_1, method = "pa")
  data_glom <- tax_glom(data_pa, taxrank = "Prevalence_alt")
  
  # Asegurarse de que los nombres de los taxones sean únicos
  taxa_names(data_glom) <- make.unique(as.character(tax_table(data_glom)[, "Prevalence_alt"]))
  
  ps_data1 <- psmelt(data_glom)
  ps_data1 <- ps_data1[, c("OTU", "Sample", "Abundance", "Prevalence_alt")]
  names(ps_data1) <- c("ASV", "Sample", "Abundance", "Prev_cat")
  ps_data1$Study <- uniq[i]
  prevalence_barplot_df[[i]] <- ps_data1
  
  # Crear la figura para cada estudio
  fig <- ps_data1 %>%
    ggplot(aes(x = Sample, y = Abundance, fill = Prev_cat)) +
    geom_bar(stat = "identity", width = 1, position = "fill") +
    scale_fill_viridis(discrete = TRUE, option = "D") +
    ggtitle(uniq[i])
  
  fig_list[[i]] <- fig
}

# Combinar las figuras individuales en una disposición 3x3
ggarrange(
  plotlist = fig_list,
  ncol = 3,
  nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)

# Unir todas las figuras en un dataframe para el gráfico combinado final
prevalence_barplot_df <- do.call(rbind, prevalence_barplot_df)
prevalence_barplot_df$Study <- factor(prevalence_barplot_df$Study, levels = uniq)
prevalence_barplot_df <- subset(prevalence_barplot_df, Prev_cat != "0")  # Eliminar taxones con prevalencia cero

# Convertir Prev_cat a factor con niveles consistentes
levels_prevalence <- c("Abajo_30%", "30-50%", "50-70%", "70-90%", "Arriba_90%")
prevalence_barplot_df$Prev_cat <- factor(prevalence_barplot_df$Prev_cat, levels = levels_prevalence)

# Crear el gráfico final combinado con etiquetas de porcentaje
G <- ggplot(prevalence_barplot_df, aes(fill = Prev_cat, y = Abundance, x = Study)) +
  geom_bar(position = "fill", stat = "identity") +
  ylab("Prevalencia promedio de ASVs") +
  theme_grey(base_size = 12) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  guides(fill = guide_legend(title = "Prevalencia")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("B") + 
  theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 16, face = "bold"))
G
################################################################################
####
# Figura 7: Curva de Abundancia-Ocupación 
H <- ggplot(occupancy_abundance_df, 
            aes(y = RelPrev, x = RelAbundance, group = host_study, col = host_study)) +
  geom_smooth(size = 2, se = F) + 
  theme_bw(base_size = 14) +
  ylim(0, 1) +
  scale_color_manual(values = palette) +
  xlab("Abundancia Relativa") +
  ylab("Prevalencia") +
  theme(legend.position = "right") +  # Incluir leyenda
  ggtitle("") +
  guides(colour = guide_legend(title = "Estudio")) +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 16, face = "bold"), legend.position = "none")
# Figura 8: Curva de Rango-Abundancia 
occupancy_abundance_df <- occupancy_abundance_df %>%
  arrange(host_study, -RelAbundance) %>%
  group_by(host_study) %>%
  mutate(Rank = rank(-RelAbundance))
I <- ggplot(occupancy_abundance_df, aes(y = RelAbundance, x = Rank, group = host_study, col = host_study)) +
  geom_smooth(size = 2, se = F) +
  theme_bw(base_size = 14) +
  xlim(0, 50) + 
  scale_color_manual(values = palette) +
  ylab("Abundancia relativa") +
  xlab("Rango de abundancia") +
  theme(legend.position = "right") +  # Incluir leyenda
  ggtitle("") +
  guides(colour = guide_legend(title = "Estudio")) + ggtitle("B") + 
  theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, size = 16, face = "bold"), legend.position = "right")
#################################
A 
#
(B + C + D) 
#
(E + G) 
#
(H + I)
##################################
# Función para filtrar ASVs que tienen asignación completa hasta el nivel especificado (ej., "")
filter_physeq_complete_taxonomy <- function(physeq_obj, min_level = "Genus") {
  # Convertir la tabla de taxonomía en un data frame
  tax_table_df <- as.data.frame(tax_table(physeq_obj))
  
  # Identificar las columnas hasta el nivel especificado
  levels <- colnames(tax_table_df)
  max_level_index <- match(min_level, levels)
  
  # Verificar que el nivel mínimo existe en la tabla de taxonomía
  if (is.na(max_level_index)) {
    stop("El nivel especificado no se encuentra en la tabla de taxonomía.")
  }
  
  # Filtrar ASVs que tienen información completa hasta el nivel mínimo especificado
  valid_asvs <- rownames(tax_table_df[complete.cases(tax_table_df[, 1:max_level_index]), ])
  
  # Prune (eliminar) ASVs que no cumplen con el criterio
  physeq_filtered <- prune_taxa(valid_asvs, physeq_obj)
  
  return(physeq_filtered)
}

# Aplicar el filtrado a cada objeto en la lista de phyloseq
phylo_list_filtered <- lapply(phylo_list, filter_physeq_complete_taxonomy, min_level = "Genus")
##################################
# Cálculo de diversidad alfa
# Lista para almacenar los resultados
list1 <- list()
uniq <- names(phylo_list_filtered)
thresholds <- c(0, 0.3, 0.5, 0.7, 0.9)

# Bucle para cada estudio
for (j in 1:length(uniq)) {
  data <- phylo_list_filtered[[j]]
  data <- prune_taxa(taxa_sums(data) > 0, data)
  
  # Corrección del argumento 'rngseed'
  data <- rarefy_even_depth(data, sample.size = 15000, rngseed = 100, replace = TRUE, trimOTUs = TRUE, verbose = FALSE)
  list2 <- list()
  
  # Bucle para cada umbral de prevalencia
  for (i in 1:length(thresholds)) {
    tryCatch({
      data_subset <- microbiome::core(data, detection = 0, prevalence = thresholds[i])
      alpha <- estimate_richness(data_subset, measures = c("Observed", "Shannon", "Chao1", "Simpson"))
      alpha$Sample <- sample_names(data_subset)
      alpha$Study <- as.character(uniq[[j]])
      alpha$Prevalence <- thresholds[i]
      alpha$SampleSums <- sample_sums(data_subset)
      list2[[i]] <- alpha
    }, error = function(e) {
      warning(paste("Error en el umbral de prevalencia", thresholds[i], "para el estudio", uniq[[j]], ":", e$message))
    })
  }
  
  # Combinar resultados del segundo bucle
  alpha_df <- do.call(rbind, list2)
  list1[[j]] <- alpha_df
}

# Combinar todos los estudios
alpha_df_all <- do.call(rbind, list1)

# Verificar que la columna 'Study' está presente antes de estandarizar
if ("Study" %in% colnames(alpha_df_all)) {
  # Crear columnas escaladas por estudio
  alpha_df_all$Observed_scaled <- ave(alpha_df_all$Observed, alpha_df_all$Study, FUN = scale)
  alpha_df_all$Chao1_scaled <- ave(alpha_df_all$Chao1, alpha_df_all$Study, FUN = scale)
  alpha_df_all$Shannon_scaled <- ave(alpha_df_all$Shannon, alpha_df_all$Study, FUN = scale)
  alpha_df_all$Simpson_scaled <- ave(alpha_df_all$Simpson, alpha_df_all$Study, FUN = scale)
}

# Asegurarte de que las posiciones de columnas son correctas antes de reordenar
# Solo seleccionar las columnas que realmente existen
existing_columns <- colnames(alpha_df_all)

# Columnas deseadas
desired_columns <- c("Sample", "Study", "Prevalence", "Observed", "Shannon", "Chao1", "Simpson", 
                     "SampleSums", "Observed_scaled", "Chao1_scaled", "Shannon_scaled", "Simpson_scaled")

# Filtrar las columnas que existen
columns_to_select <- intersect(desired_columns, existing_columns)

# Seleccionar solo las columnas existentes
alpha_df_all <- alpha_df_all[, columns_to_select]

# Revisar la estructura final
str(alpha_df_all)
#####################
alpha_df_all<-alpha_df_all[,c(5,6,7, 1,2,3,4,9,10,11, 12, 8)]
#Vizualizacion
# Preparación de los datos para diversidad alfa escalada
alpha_df_all$Study <- factor(alpha_df_all$Study, levels = uniq)
alpha_df_all$Prevalence <- factor(alpha_df_all$Prevalence)
# Seleccionar las columnas necesarias para los cálculos de diversidad alfa escalada
alpha_short_scaled <- alpha_df_all[, c("Sample", "Study", "Prevalence", "Observed", "Chao1", "Shannon", "Simpson")]
# Convertir los datos a formato largo para incluir todas las métricas
alpha_long_scaled <- gather(alpha_short_scaled, Index, Distance, Observed, Chao1, Shannon, Simpson, factor_key = TRUE)
# Asegurar que los niveles de "Index" están correctamente especificados
alpha_long_scaled$Index <- factor(alpha_long_scaled$Index, levels = c("Observed", "Chao1", "Shannon", "Simpson"))
# Visualización de la diversidad alfa escalada
J <- ggplot(alpha_long_scaled, aes(x = Prevalence, y = Distance, fill = Study, col = Study, group = interaction(Study, Index))) +
  geom_jitter(pch = 21, size = 0.5, alpha = 0.7, width = 0.3) +
  facet_wrap(~ Study + Index, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 14) +
  xlab("Prevalencia") +
  ylab("Alfa diversidad") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  ggtitle("Diversidad Alfa Escalada")
# Preparación de los datos para diversidad alfa no escalada
alpha_short_unscaled <- alpha_df_all[, c("Sample", "Study", "Prevalence", "Observed", "Chao1", "Shannon", "Simpson")]
# Convertir los datos a formato largo para incluir todas las métricas no escaladas
alpha_long_unscaled <- gather(alpha_short_unscaled, Index, Distance, Observed, Chao1, Shannon, Simpson, factor_key = TRUE)
# Asegurar que los niveles de "Index" están correctamente especificados
alpha_long_unscaled$Index <- factor(alpha_long_unscaled$Index, levels = c("Observed", "Chao1", "Shannon", "Simpson"))
# Visualización de la diversidad alfa no escalada
K <- ggplot(alpha_long_unscaled, aes(x = Prevalence, y = Distance, fill = Study, col = Study, group = interaction(Study, Index))) +
  geom_jitter(pch = 21, size = 0.5, alpha = 0.7, width = 0.3) +
  facet_wrap(~ Study + Index, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 14) +
  xlab("Prevalencia") +
  ylab("Alfa Diversidad") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  ggtitle("Diversidad Alfa Sin Escalar")
####################
J
#
K
#####
######
####
# Calculo de diversidad Beta 
thresholds <- c(0, 0.3, 0.5, 0.7, 0.9)
# Lista para almacenar resultados de PCoA por prevalencia
pcoa_results_list <- list()
# Lista para almacenar data frames de distancias por prevalencia
beta_df_by_prevalence <- list()

# Iterar sobre cada umbral de prevalencia
for (i in 1:length(thresholds)) {
  prev_threshold <- thresholds[i]
  
  list1 <- list()
  uniq <- names(phylo_list_filtered)
  
  # For loop para cada estudio
  for (j in 1:length(uniq)) {
    data <- phylo_list_filtered[[j]]  # for loop 1 (uniq)
    data <- prune_taxa(taxa_sums(data) > 0, data)
    data <- rarefy_even_depth(data, sample.size = 15000, rngsee = 100, replace = TRUE,
                              trimOTUs = TRUE, verbose = TRUE)
    
    tryCatch({ # catch errors
      # Filtrar datos por umbral de prevalencia
      data_subset <- microbiome::core(data, detection = 0, prevalence = prev_threshold)
      
      # Calcular distancias
      unifrac <- as.matrix(phyloseq::distance(data_subset, method = "unifrac"))
      wunifrac <- as.matrix(phyloseq::distance(data_subset, method = "wunifrac"))
      morisita <- as.matrix(phyloseq::distance(data_subset, method = "morisita"))
      jaccard <- as.matrix(phyloseq::distance(data_subset, method = "jaccard"))
      braycurtis <- as.matrix(phyloseq::distance(data_subset, method = "bray"))  
      
      # Crear el dataframe con las distancias
      beta <- data.frame(colMeans(unifrac, na.rm = TRUE))
      names(beta) <- "unifrac"
      beta$wunifrac <- colMeans(wunifrac, na.rm = TRUE)
      beta$morisita <- colMeans(morisita, na.rm = TRUE)
      beta$jaccard <- colMeans(jaccard, na.rm = TRUE)
      beta$braycurtis <- colMeans(braycurtis, na.rm = TRUE)  # Bray-Curtis
      beta$Sample <- sample_names(data_subset)
      beta$Study <- as.character(uniq[[j]])
      beta$Prevalence <- prev_threshold
      beta$SampleSums <- sample_sums(data_subset)
      beta$ASV_Abundance <- rowSums(otu_table(data_subset))
      beta$Num_ASVs <- rowSums(otu_table(data_subset) > 0)
      
      list1[[j]] <- beta
    }, error = function(e) {})
  }
  
  # Combinar todos los estudios en un solo data frame para este umbral de prevalencia
  beta_df <- do.call(rbind, list1)
  
  # Guardar el data frame en la lista por prevalencia
  beta_df_by_prevalence[[i]] <- beta_df
  
  # Realizar el análisis de PCoA (o PSA) para cada umbral
  distance_matrix <- as.dist(beta_df$braycurtis)  # Cambia 'braycurtis' por la distancia que prefieras
  
  # Análisis de PCoA
  pcoa_result <- pcoa(distance_matrix)
  
  # Extraer los ejes principales
  ordination_df <- data.frame(pcoa_result$vectors[, 1:2])  # Tomar los primeros dos ejes
  colnames(ordination_df) <- c("Axis.1", "Axis.2")  # Renombrar las columnas para evitar confusiones
  ordination_df$Study <- beta_df$Study
  ordination_df$Prevalence <- as.factor(beta_df$Prevalence)
  
  # Guardar el resultado del PCoA en la lista
  pcoa_results_list[[i]] <- ordination_df
  
  # Asegúrate de que la columna ASV_Abundance esté correctamente definida y pertenece al mismo dataframe que las distancias
  beta_df$ASV_Abundance <- as.numeric(beta_df$Num_ASVs)
  
  # Realizar el PERMANOVA asegurándote de que el dataframe beta_df esté siendo utilizado en el análisis
  permanova_result <- adonis2(distance_matrix ~ Num_ASVs, data = beta_df, permutations = 999)
  
  # Extraer el valor de p y R² del resultado de PERMANOVA
  r_squared <- permanova_result$R2[1]
  p_value <- permanova_result$`Pr(>F)`[1]
  
  # Definir el número de asteriscos según el valor de p
  sig_level <- ifelse(p_value <= 0.001, "***",
                      ifelse(p_value <= 0.01, "**",
                             ifelse(p_value <= 0.05, "*", "")))
  
  # Modificar el gráfico para incluir los resultados de PERMANOVA
 L <- ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, color = Study)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(geom = "path", aes(group = Study), size = 1, linetype = "solid") +
    scale_color_manual(values = palette) +
    labs(
      x = paste0("PCoA Axis 1 (", round(pcoa_result$values$Relative_eig[1] * 100, 2), "%)"),
      y = paste0("PCoA Axis 2 (", round(pcoa_result$values$Relative_eig[2] * 100, 2), "%)"),
      color = "Study"
    ) +
    theme_minimal(base_size = 14) +
    ggtitle(paste0("PCoA de diversidad beta: Prevalencia ", thresholds[i], 
                   "\nR² = ", round(r_squared, 4), ", p = ", round(p_value, 4)," ", sig_level)) +
    theme(legend.position = "right") +
    guides(colour = guide_legend(title = "Estudio"))
  
  # Imprimir el gráfico para visualizarlo
  print(L)
}
################################################################################
# Diversidad Beta para comparar metricas de distancia
# Lista para almacenar los gráficos por prevalencia
plots_by_prevalence <- list()

# Iterar sobre cada umbral de prevalencia
for (i in 1:length(thresholds)) {
  prev_threshold <- thresholds[i]
  
  list1 <- list()
  uniq <- names(phylo_list_filtered)
  
  for (j in 1:length(uniq)) {
    data <- phylo_list_filtered[[j]]
    data <- prune_taxa(taxa_sums(data) > 0, data)
    data <- rarefy_even_depth(data, sample.size = 15000, rngsee = 100, replace = TRUE,
                              trimOTUs = TRUE, verbose = TRUE)
    
    tryCatch({
      data_subset <- microbiome::core(data, detection = 0, prevalence = prev_threshold)
      
      # Calcular distancias
      braycurtis <- as.matrix(phyloseq::distance(data_subset, method = "bray"))
      morisita <- as.matrix(phyloseq::distance(data_subset, method = "morisita"))
      
      # Crear data frames para Bray-Curtis y Morisita
      beta_bray <- data.frame(colMeans(braycurtis, na.rm = TRUE))
      names(beta_bray) <- "braycurtis"
      beta_bray$Study <- as.character(uniq[[j]])
      beta_bray$Prevalence <- prev_threshold
      beta_bray$Sample <- sample_names(data_subset)
      beta_bray$SampleSums <- sample_sums(data_subset)
      beta_bray$ASV_Abundance <- rowSums(otu_table(data_subset))
      beta_bray$Num_ASVs <- rowSums(otu_table(data_subset) > 0)
      beta_morisita <- data.frame(colMeans(morisita, na.rm = TRUE))
      names(beta_morisita) <- "morisita"
      beta_morisita$Study <- as.character(uniq[[j]])
      beta_morisita$Prevalence <- prev_threshold
      beta_morisita$Sample <- sample_names(data_subset)
      beta_morisita$SampleSums <- sample_sums(data_subset)
      beta_morisita$ASV_Abundance <- rowSums(otu_table(data_subset))
      beta_morisita$Num_ASVs <- rowSums(otu_table(data_subset) > 0)
      
      list1[[j]] <- list(bray = beta_bray, morisita = beta_morisita)
    }, error = function(e) {})
  }
  
  # Combinar data frames de todos los estudios
  beta_bray_df <- do.call(rbind, lapply(list1, function(x) x$bray))
  beta_morisita_df <- do.call(rbind, lapply(list1, function(x) x$morisita))
  
  # Crear matrices de distancias para PCoA
  distance_matrix_bray <- as.dist(beta_bray_df$braycurtis)
  distance_matrix_morisita <- as.dist(beta_morisita_df$morisita)
  
  # Realizar PCoA
  pcoa_bray <- pcoa(distance_matrix_bray)
  pcoa_morisita <- pcoa(distance_matrix_morisita)
  
  # Data frames para PCoA
  ordination_bray <- data.frame(pcoa_bray$vectors[, 1:2])
  colnames(ordination_bray) <- c("Axis.1", "Axis.2")
  ordination_bray$Study <- beta_bray_df$Study
  ordination_bray$Prevalence <- as.factor(beta_bray_df$Prevalence)
  ordination_bray$Distance <- "Bray-Curtis"
  
  ordination_morisita <- data.frame(pcoa_morisita$vectors[, 1:2])
  colnames(ordination_morisita) <- c("Axis.1", "Axis.2")
  ordination_morisita$Study <- beta_morisita_df$Study
  ordination_morisita$Prevalence <- as.factor(beta_morisita_df$Prevalence)
  ordination_morisita$Distance <- "Morisita"
  
  # Combinar data frames para gráficos
  ordination_combined <- rbind(ordination_bray, ordination_morisita)
  
  # Asegúrate de que la columna ASV_Abundance esté correctamente definida y pertenece al mismo dataframe que las distancias
  beta_bray_df$ASV_Abundance <- as.numeric(beta_bray_df$Num_ASVs)
  beta_morisita_df$ASV_Abundance <- as.numeric(beta_morisita_df$Num_ASVs)
  
  # Resultados de PERMANOVA
  permanova_bray <- adonis2(distance_matrix_bray ~ Num_ASVs, data = beta_bray_df, permutations = 999)
  permanova_morisita <- adonis2(distance_matrix_morisita ~ Num_ASVs, data = beta_morisita_df, permutations = 999)
  
  r_squared_bray <- permanova_bray$R2[1]
  p_value_bray <- permanova_bray$`Pr(>F)`[1]
  r_squared_morisita <- permanova_morisita$R2[1]
  p_value_morisita <- permanova_morisita$`Pr(>F)`[1]
  
  sig_level_bray <- ifelse(p_value_bray <= 0.001, "***",
                           ifelse(p_value_bray <= 0.01, "**",
                                  ifelse(p_value_bray <= 0.05, "*", "")))
  
  sig_level_morisita <- ifelse(p_value_morisita <= 0.001, "***",
                               ifelse(p_value_morisita <= 0.01, "**",
                                      ifelse(p_value_morisita <= 0.05, "*", "")))
  
  # Crear gráficos de Bray-Curtis y Morisita
  plot_bray <- ggplot(subset(ordination_combined, Distance == "Bray-Curtis"), 
                      aes(x = Axis.1, y = Axis.2, color = Study)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(geom = "path", aes(group = Study), size = 1, linetype = "solid") +
    labs(
      x = paste0("PCoA Axis 1 (", round(pcoa_bray$values$Relative_eig[1] * 100, 2), "%)"),
      y = paste0("PCoA Axis 2 (", round(pcoa_bray$values$Relative_eig[2] * 100, 2), "%)"),
      title = paste0("PCoA: Prevalencia ", thresholds[i], "\nBray-Curtis: R² = ", 
                     round(r_squared_bray, 4), ", p = ", round(p_value_bray, 4), " ", sig_level_bray),
      color = "Estudio"
    ) +
    theme_minimal(base_size = 14) +
    scale_color_manual(values = palette)
  
  plot_morisita <- ggplot(subset(ordination_combined, Distance == "Morisita"), 
                          aes(x = Axis.1, y = Axis.2, color = Study)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(geom = "path", aes(group = Study), size = 1, linetype = "solid") +
    labs(
      x = paste0("PCoA Axis 1 (", round(pcoa_morisita$values$Relative_eig[1] * 100, 2), "%)"),
      y = paste0("PCoA Axis 2 (", round(pcoa_morisita$values$Relative_eig[2] * 100, 2), "%)"),
      title = paste0("\nMorisita: R² = ", 
                     round(r_squared_morisita, 4), ", p = ", round(p_value_morisita, 4), " ", sig_level_morisita),
      color = "Estudio"
    ) +
    theme_minimal(base_size = 14) +
    scale_color_manual(values = palette)
  
  # Combinar gráficos en un solo panel
  combined_plot <- plot_bray + plot_morisita + plot_layout(guides = "collect") & 
    theme(legend.position = "right")
  
  # Guardar el gráfico en la lista
  plots_by_prevalence[[i]] <- combined_plot
  
  # Mostrar el gráfico
  print(combined_plot)
}
#################################################################################################
# Vizualizacion de Abundancia
# Definir los umbrales de prevalencia con los que trabajamos
prevalence_thresholds <- c(0, 0.3, 0.5, 0.7, 0.9)
# Calculo de abundancias 
get_abundances_by_threshold <- function(phylo_list_filtered, threshold_list, taxonomic_rank, detection = 0.001) {
  abundances_list <- list()
  
  for (i in 1:length(phylo_list_filtered)) {
    physeq_obj <- phylo_list_filtered[[i]]
    
    list_abundances <- list()
    
    for (j in 1:length(threshold_list)) {
      # Estimar el core microbioma según el umbral de prevalencia
      data_core <- microbiome::core(physeq_obj, detection = detection, prevalence = threshold_list[j])
      
      # Extraer tabla de taxonomía y abundancias
      tax_table_core <- tax_table(physeq_obj)[taxa_names(physeq_obj) %in% taxa_names(data_core), ]
      
      # Comprobar si el nivel taxonómico existe en la tabla
      if (!taxonomic_rank %in% colnames(tax_table_core)) {
        warning(paste("El nivel taxonómico", taxonomic_rank, "no se encuentra en la tabla de taxonomía"))
        next
      }
      
      # Calcular abundancias de taxones
      taxon_abundances <- taxa_sums(data_core)
      taxon_names <- as.character(tax_table_core[, taxonomic_rank])
      
      # Si hay NAs, reemplazarlos por "Other" para evitar problemas en la visualización
      taxon_names[is.na(taxon_names)] <- "Other"
      
      # Crear un data frame con abundancias
      abundances_df <- data.frame(
        Taxon = taxon_names,
        Abundance = taxon_abundances,
        Prevalence = threshold_list[j],
        Study = names(phylo_list_filtered)[i]
      )
      
      list_abundances[[j]] <- abundances_df
    }
    
    abundances_list[[i]] <- do.call(rbind, list_abundances)
  }
  
  return(do.call(rbind, abundances_list))
}
# Obtener abundancias a nivel de filo (Phylum)
abundances_phylum <- get_abundances_by_threshold(phylo_list_filtered, prevalence_thresholds, "Phylum", detection = 0.001)
# Definir una paleta personalizada de 31 colores
custom_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF0F", "#984EA3", "#FF7F00", 
  "#FFFF33", "#A65628", "#F781BF", "darkgreen", "#66C2A5",
  "#FC8D30", "#8DA0CB", "#E78A50", "#A6D854", "#FFD92F",
  "#E5C450", "darkred", "#8DD3C7", "#FFFFC9", "#BEBADA", 
  "#FB8072", "#80B100", "#FDB462", "darkblue", "#FCCDE5",
  "#D9D970", "#BC80BD", "#CCEBC5", "#FFED6F", "#FFB3B3", 
  "#Baa"
)
# Crear el gráfico con la paleta personalizada
M <- ggplot(abundances_phylum, aes(x = factor(Prevalence), y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_palette) +  # Aplicar la paleta manual
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Umbral de prevalencia", y = "Número de ASVs", fill = "Phylum") +
  ylim(0, 50000000) +
  theme_minimal(base_size = 14) +
  ggtitle("N° de ASVs por prevalencia (Abundancia a nivel de filo)") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.key.size = unit(0.4, 'cm'),
    axis.line = element_line(color = "black"),  # Añadir líneas de los ejes
    axis.ticks = element_line(color = "black"),  # Añadir marcadores en los ejes
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(title = "Filo", ncol = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # Ajustar la leyenda a dos columnas

M
#####################################################
# Vizualizacion de abundancia por separado 
# Definir los límites del eje Y manualmente para cada umbral
ylim_values <- list(
  "0"   = 50000000,  # Para prevalencia 0
  "0.3" = 35000000,  # Ajuste manual para prevalencia 0.3
  "0.5" = 25000000,   # Ajuste manual para prevalencia 0.5
  "0.7" = 16000000,   # Ajuste manual para prevalencia 0.7
  "0.9" = 7000000    # Ajuste manual para prevalencia 0.9
)

# Inicializar lista para almacenar gráficos por prevalencia
plots_by_prevalence <- list()

# Paso 2: Crear gráficos por cada umbral de prevalencia y ajustar el eje Y manualmente
for (prev in prevalence_thresholds) {
  # Filtrar el dataframe para el umbral actual
  df_prev <- abundances_phylum[abundances_phylum$Prevalence == prev, ]
  
  # Comprobar si hay valores de abundancia antes de continuar
  if (nrow(df_prev) == 0) {
    warning(paste("No data available for prevalence:", prev))
    next
  }
  
  # Obtener el valor del eje Y para este umbral
  max_y_value <- ylim_values[[as.character(prev)]]
  
  # Crear el gráfico ajustando el eje Y manualmente
  plot <- ggplot(df_prev, aes(x = factor(Prevalence), y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = custom_palette) +  # Aplicar la paleta manual
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Umbral de prevalencia", y = "Número de ASVs", fill = "Phylum") +
    ylim(0, max_y_value) +  # Ajuste manual del límite del eje Y
    theme_minimal(base_size = 14) +
    ggtitle(paste0("N° de ASVs (Abundancia a nivel de filo) - Prevalencia ", prev)) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      legend.key.size = unit(0.4, 'cm'),
      axis.line = element_line(color = "black"),  # Añadir líneas de los ejes
      axis.ticks = element_line(color = "black"),  # Añadir marcadores en los ejes
      legend.text = element_text(size = 12),
      legend.position = "right"
    ) +
    guides(fill = guide_legend(title = "Filo", ncol = 1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # Ajustar la leyenda a dos columnas
  
  # Almacenar el gráfico en la lista
  plots_by_prevalence[[as.character(prev)]] <- plot
}

# Paso 3: Visualizar cada gráfico por prevalencia
for (plot in plots_by_prevalence) {
  print(plot)
}

################################################################################
# Core Microbiome 
# Visualizacion taxonomica por muestra (Estudio)
# Definir listas para almacenar los resultados de cada análisis
p_heatmap_asvs_list <- list()
p_heatmap_genus_list <- list()
p_lineplot_list <- list()
# Definir los umbrales de prevalencia con los que trabajamos
prevalence_thresholds <- c(0, 0.3, 0.5, 0.7, 0.9)

# Iterar a través de cada objeto phyloseq en phylo_list
for (i in 1:length(phylo_list_filtered)) {
  physeq_obj <- phylo_list_filtered[[i]]
  
  # Mantener solo taxones con sumas positivas
  physeq_obj <- prune_taxa(taxa_sums(physeq_obj) > 0, physeq_obj)
  
  # Calcular la versión compositiva de los datos (abundancias relativas)
  physeq_obj_compositional <- microbiome::transform(physeq_obj, "compositional")
  
  # Asignar secuencias ASVs si están en formato de secuencia de nucleótidos
  physeq_obj_compositional <- microbiome::add_refseq(physeq_obj_compositional)
  
  # Limpiar nombres de taxones
  tax_table(physeq_obj_compositional)[, colnames(tax_table(physeq_obj_compositional))] <- 
    gsub(tax_table(physeq_obj_compositional)[, colnames(tax_table(physeq_obj_compositional))], 
         pattern = "[a-z]__", replacement = "")
  
  # Line plot de core microbiome a nivel de ASVs
  detections <- c(0.0001, 0.0005, 0.001, 0.01, 0.1, 1, 10)
  prevalence_t <- seq(.05, 1, .05)
  # Ajuste de detección menos conservador
  p_lineplot <- plot_core(physeq_obj_compositional, prevalences = prevalence_t, detections = detections, plot.type = "lineplot") + 
    xlab("Abundancia relativa (%)") + 
    ylab("Tamaño del núcleo central (N)") +
    guides(color = guide_legend(title = "Prevalencia")) +
    theme_bw() + 
    ggtitle(paste("Grafico de lineas centrales ASVs - Estudio", i))
  
  # Guardar el gráfico de line plot
  p_lineplot_list[[i]] <- p_lineplot
  
  # Heatmap core microbiome a nivel de ASVs
  detections <- round(10^seq(log10(1e-4), log10(.1), length = 10), 3)  # Ajustar para incluir más ASVs
  p_heatmap_asvs <- tryCatch({
    plot_core(physeq_obj_compositional, 
              plot.type = "heatmap", 
              prevalences = prevalence_thresholds, 
              detections = detections, min.prevalence = .3,
              colours = sequential_hcl(n = 10, palette = "Inferno", rev = TRUE)) +
      xlab("Umbral de detección (abundancia relativa (%))") + 
      ylab("ASVs") +
      guides(fill = guide_legend(title = "Prevalencia", reverse = TRUE)) +  # Leyenda invertida
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1),
            axis.text.y = element_text(face = "italic")) + 
      ggtitle(paste("Heatmap ASVs - Estudio", i))
  }, error = function(e) {
    warning(paste("Error al generar el heatmap de ASVs para el estudio", i, ":", e$message))
    NULL
  })
  
  # Guardar el heatmap de ASVs si se generó correctamente
  if (!is.null(p_heatmap_asvs)) {
    p_heatmap_asvs_list[[i]] <- p_heatmap_asvs
  }
  
  # A nivel de género
  physeq_obj_genus <- aggregate_taxa(physeq_obj_compositional, "Genus") # Agrupar al nivel de género
  
  # Eliminar taxones clasificados como "Unknown"
  physeq_obj_genus <- subset_taxa(physeq_obj_genus, Genus != "Unknown")
  
  # Heatmap del core microbiome a nivel de género
  p_heatmap_genus <- tryCatch({
    plot_core(physeq_obj_genus, 
              plot.type = "heatmap",
              prevalences = prevalence_thresholds, 
              detections = detections, min.prevalence = .3, 
              colours = sequential_hcl(n = 10, palette = "Inferno", rev = TRUE)) +
      xlab("Umbral de detección (abundancia relativa (%))") + 
      ylab("Géneros") + 
      guides(fill = guide_legend(title = "Prevalencia", reverse = TRUE)) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1),
            axis.text.y = element_text(face = "italic")) + 
      ggtitle(paste("Heatmap Género - Estudio", i)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }, error = function(e) {
    warning(paste("Error al generar el heatmap de géneros para el estudio", i, ":", e$message))
    NULL
  })
  
  # Guardar el heatmap de géneros si se generó correctamente
  if (!is.null(p_heatmap_genus)) {
    p_heatmap_genus_list[[i]] <- p_heatmap_genus
  }
}

# Filtrar gráficos válidos en p_lineplot_list, p_heatmap_asvs_list y p_heatmap_genus_list
p_lineplot_list <- Filter(Negate(is.null), p_lineplot_list)
p_heatmap_asvs_list <- Filter(Negate(is.null), p_heatmap_asvs_list)
p_heatmap_genus_list <- Filter(Negate(is.null), p_heatmap_genus_list)

# Visualizar los gráficos de line plot para cada estudio
for (i in 1:length(p_lineplot_list)) {
  print(p_lineplot_list[[i]])
}

# Visualizar los heatmaps de ASVs para cada estudio
for (i in 1:length(p_heatmap_asvs_list)) {
  print(p_heatmap_asvs_list[[i]])
}

# Visualizar los heatmaps de géneros para cada estudio
for (i in 1:length(p_heatmap_genus_list)) {
  print(p_heatmap_genus_list[[i]])
}

################################################################################
# Visualizacion Taxonomica por Prevalencia 
# Función para combinar y alinear géneros entre estudios
combine_abundances <- function(abundance_list) {
  combined_df <- Reduce(function(x, y) {
    full_join(x, y, by = "Taxon")
  }, abundance_list)
  
  # Reemplazar los NA con 0
  combined_df[is.na(combined_df)] <- 0
  rownames(combined_df) <- combined_df$Taxon
  combined_df$Taxon <- NULL  # Remover columna Taxon
  return(combined_df)
}

# Filtrar los géneros más prevalentes y ajustar el heatmap
for (prev_threshold in prevalence_thresholds) {
  
  genera_list <- list()
  
  for (i in 1:length(phylo_list_filtered)) {
    physeq_obj <- phylo_list_filtered[[i]]
    
    # Transformación composicional
    physeq_obj_compositional <- microbiome::transform(physeq_obj, "compositional")
    
    # Agrupar a nivel de género
    physeq_obj_genus <- aggregate_taxa(physeq_obj_compositional, "Phylum") # Configurar nivel taxonomico segun prefiera 
    
    # Eliminar taxones clasificados como "Unknown"
    physeq_obj_genus <- subset_taxa(physeq_obj_genus, Phylum != "Unknown") # Configurar nivel taxonomico segun prefiera 
    
    # Filtrar géneros según el umbral de prevalencia
    core_genera <- core(physeq_obj_genus, detection = 0.0001, prevalence = prev_threshold)
    
    # Extraer tabla de abundancia para los géneros seleccionados
    abundances <- abundances(core_genera)
    
    # Convertir la tabla a un data.frame e incluir los nombres de los taxones
    df_abundances <- data.frame(Taxon = rownames(abundances), abundances)
    genera_list[[i]] <- df_abundances
  }
  
  # Combinar los resultados de todos los estudios alineando los géneros
  combined_abundances <- combine_abundances(genera_list)
  
  # Filtrar géneros con abundancia acumulada baja (opcional para simplificar el gráfico)
  total_abundance <- rowSums(combined_abundances)
  combined_abundances <- combined_abundances[total_abundance > quantile(total_abundance, 0.30), ]  # Filtrar por el 25% más abundante
  
  # Generar el heatmap ajustando la escala y el texto
  pheatmap(combined_abundances, 
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           main = paste("Heatmap - Prevalencia:", prev_threshold),
           scale = "none",  # Escalar por fila
           color = sequential_hcl(n = 10, palette = "Inferno", rev = TRUE),
           fontsize_row = 9,  # Tamaño de fuente para las filas
           fontsize_col = 0.0001,  # Tamaño de fuente para las columna, se elige un numero bajo para no saturar el plot con las etiquetas de las muestras 
           show_colnames = TRUE,  # Ocultarlo solo como opcion 
           angle_col = 45)  # Rotar etiquetas de columna para mejorar visibilidad
}
################################################################################
# Core Microbiome - Heatmap Final - Nivel de Filo
# Parámetros de detección y prevalencia
detections <- seq(0.0002, 0.1, by = 0.005)
prevalences <- c(0, 0.3, 0.5, 0.7, 0.9)
# Función modificada para combinar datos de todos los estudios y generar un solo heatmap
plot_combined_heatmap_by_prevalence <- function(phylo_list_filtered, taxonomic_level = "Phylum") {
  
  combined_physeq <- NULL  # Inicializar variable para almacenar los datos combinados
  
  # Iterar sobre los objetos phyloseq en la lista
  for (i in 1:length(phylo_list_filtered)) {
    physeq_obj <- phylo_list_filtered[[i]]
    
    # Convertir a abundancias relativas
    physeq_obj_compositional <- microbiome::transform(physeq_obj, "compositional")
    
    # Agrupar a nivel taxonómico
    physeq_obj_phylum <- aggregate_taxa(physeq_obj_compositional, taxonomic_level)
    
    # Eliminar taxones clasificados como "Unknown"
    physeq_obj_phylum <- subset_taxa(physeq_obj_phylum, Phylum != "Unknown")
    
    # Combinar los objetos phyloseq (si hay más de uno)
    if (is.null(combined_physeq)) {
      combined_physeq <- physeq_obj_phylum
    } else {
      combined_physeq <- merge_phyloseq(combined_physeq, physeq_obj_phylum)
    }
  }
  
  # Generar el heatmap con los datos combinados
  p <- plot_core(combined_physeq, 
                 plot.type = "heatmap", 
                 colours = sequential_hcl(n = 10, palette = "Rocket", rev = TRUE),  # Colores para el heatmap
                 prevalences = prevalences, 
                 detections = detections, 
                 min.prevalence = 0.3) +  # Cambia este valor según lo que quieras ver
    xlab("Umbral de detección (abundancia relativa (%))") +
    ylab(taxonomic_level) +
    labs(y = "Filo") +
    guides(fill = guide_legend(title = "Prevalencia", reverse = TRUE)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text.y = element_text(face = "italic")) +
    ggtitle(paste("Núcleo Central - Nivel de Filo")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Mostrar el heatmap combinado
  print(p)
}

# Ejecutar la función para generar un solo heatmap
plot_combined_heatmap_by_prevalence(phylo_list_filtered, "Phylum")
#####################################################################
# Core Microbiome - Heatmap Final - Nivel de Genero 
# Parámetros de detección y prevalencia
detections <- seq(0.0002, 0.008, by = 0.0002)  # De 0 a 10% en pasos de 1%
prevalences <- c(0, 0.3, 0.5, 0.7, 0.9)
# Función modificada para combinar datos de todos los estudios y generar un solo heatmap
plot_combined_heatmap_by_prevalence <- function(phylo_list_filtered, taxonomic_level = "Genus") {
  
  combined_physeq <- NULL  # Inicializar variable para almacenar los datos combinados
  
  # Iterar sobre los objetos phyloseq en la lista
  for (i in 1:length(phylo_list_filtered)) {
    physeq_obj <- phylo_list_filtered[[i]]
    
    # Convertir a abundancias relativas
    physeq_obj_compositional <- microbiome::transform(physeq_obj, "compositional")
    
    # Agrupar a nivel taxonómico
    physeq_obj_genus <- aggregate_taxa(physeq_obj_compositional, taxonomic_level)
    
    # Eliminar taxones clasificados como "Unknown"
    physeq_obj_genus <- subset_taxa(physeq_obj_genus, Genus != "Unknown")
    
    # Combinar los objetos phyloseq (si hay más de uno)
    if (is.null(combined_physeq)) {
      combined_physeq <- physeq_obj_genus
    } else {
      combined_physeq <- merge_phyloseq(combined_physeq, physeq_obj_genus)
    }
  }
  
  # Generar el heatmap con los datos combinados
  p <- plot_core(combined_physeq, 
                 plot.type = "heatmap", 
                 colours = sequential_hcl(n = 10, palette = "Inferno", rev = TRUE),  # Colores para el heatmap
                 prevalences = prevalences, 
                 detections = detections, 
                 min.prevalence = 0.3) +  # Cambia este valor según lo que quieras ver
    xlab("Umbral de detección (abundancia relativa (%))") +
    ylab(taxonomic_level) +
    labs(y = "Género") +
    guides(fill = guide_legend(title = "Prevalencia", reverse = TRUE)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text.y = element_text(face = "italic")) +
    ggtitle(paste("Núcleo Central - Nivel de Género")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  # Mostrar el heatmap combinado
  print(p)
}

# Ejecutar la función para generar un solo heatmap
plot_combined_heatmap_by_prevalence(phylo_list_filtered, "Genus")
##########################################################
# Core Microbiome - Heatmap Final - Nivel de Especie 
# Función para limpiar los nombres de las especies, eliminando los códigos de secuencia
clean_species_names <- function(physeq_obj) {
  # Limpiar nombres de especies, eliminando todo después de "_" en el nombre
  tax_table(physeq_obj)[, "Species"] <- gsub("\\(.*\\)", "", tax_table(physeq_obj)[, "Species"]) # Eliminar texto entre paréntesis
  return(physeq_obj)
}
# Parámetros de detección y prevalencia
detections <- round(10^seq(log10(1e-4), log10(.1), length = 10), 5)
prevalences <- c(0, 0.3, 0.5, 0.7, 0.9)

# Función modificada para combinar datos de todos los estudios y generar un solo heatmap
plot_combined_heatmap_by_prevalence <- function(phylo_list_filtered, taxonomic_level = "Species") {
  
  combined_physeq <- NULL  # Inicializar variable para almacenar los datos combinados
  
  # Iterar sobre los objetos phyloseq en la lista
  for (i in 1:length(phylo_list_filtered)) {
    physeq_obj <- phylo_list_filtered[[i]]
    
    # Limpiar los nombres de las especies para eliminar los códigos
    physeq_obj <- clean_species_names(physeq_obj)
    
    # Convertir a abundancias relativas
    physeq_obj_compositional <- microbiome::transform(physeq_obj, "compositional")
    
    # Agrupar a nivel taxonómico (en este caso "Species")
    physeq_obj_species <- aggregate_taxa(physeq_obj_compositional, taxonomic_level)
    
    # Eliminar taxones clasificados como "Unknown"
    physeq_obj_species <- subset_taxa(physeq_obj_species, Species != "Unknown")
    
    # Combinar los objetos phyloseq (si hay más de uno)
    if (is.null(combined_physeq)) {
      combined_physeq <- physeq_obj_species
    } else {
      combined_physeq <- merge_phyloseq(combined_physeq, physeq_obj_species)
    }
  }
  
  # Generar el heatmap con los datos combinados
  p <- plot_core(combined_physeq, 
                 plot.type = "heatmap", 
                 colours = sequential_hcl(n = 10, palette = "Inferno", rev = TRUE), # Colores para el heatmap
                 prevalences = prevalences, 
                 detections = detections, 
                 min.prevalence = 0.3) +  # Cambia este valor según lo que quieras ver
    xlab("Umbral de detección (abundancia relativa (%))") +
    ylab(taxonomic_level) +
    labs(y = "Especie") +
    guides(fill = guide_legend(title = "Prevalencia", reverse = TRUE)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text.y = element_text(face = "italic")) +
    ggtitle(paste("Núcleo Central - Nivel de Especie")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Mostrar el heatmap combinado
  print(p)
}

# Ejecutar la función para generar un solo heatmap a nivel de especie
plot_combined_heatmap_by_prevalence(phylo_list_filtered, "Species")
#####################################################################################################