
# descriptive_analysis.R
# Script de análisis descriptivo y de redes para Urban Mobilization Diffusion
# ---------------------------------------------------------------------------
# Librerías y setup
suppressMessages({
  library(tidyverse)
  library(igraph)
  library(tidygraph)
  library(statnet)
  library(ergm)
  library(ggplot2)
  library(network)
  library(kableExtra)
  library(plotly)
  library(gridExtra)
  library(viridis)
  library(sna)
  library(concorR)
  library(ggraph)
  library(Matrix)
  library(patchwork)
  library(texreg)
})

# ---------------------------
# Funciones auxiliares
# ---------------------------

# Función para leer matrices de red (archivo CSV sin header)
read_network_matrix <- function(path) {
  if(!file.exists(path)) {
    stop(paste("El archivo", path, "no existe."))
  }
  as.matrix(read.csv(path, header = FALSE))
}

# Función para calcular estadísticas de red a partir de una matriz
network_stats <- function(mat, mode = "directed") {
  g <- graph_from_adjacency_matrix(mat, mode = mode)
  tibble(
    Densidad = edge_density(g),
    N_lazos = sum(mat),
    N_nodos_activos = sum(rowSums(mat) > 0),
    Reciprocidad = reciprocity(g),
    Transitivity = transitivity(g),
    Componentes = count_components(g),
    Diametro = diameter(g, directed=TRUE, weights=NA)
  )
}

# ---------------------------
# Lectura y preparación de datos
# ---------------------------
cat("Leyendo archivos de datos...\n")
# Definir rutas
mats_paths <- list(
  Coop = "data/CoopNet.csv",
  Confianza = "data/ConfianzaNet.csv",
  Recursos = "data/RecursosNet.csv",
  Valores = "data/ValoresNet.csv",
  Parentesco = "data/ParentescoNet.csv"
)

# Leer matrices
coop_matrix       <- read_network_matrix(mats_paths$Coop)
trust_matrix      <- read_network_matrix(mats_paths$Confianza)
resources_matrix  <- read_network_matrix(mats_paths$Recursos)
values_matrix     <- read_network_matrix(mats_paths$Valores)
parentesco_matrix <- read_network_matrix(mats_paths$Parentesco)

# Leer atributos (asegúrate que no existan espacios extra en el nombre de columnas)
attributes <- read.csv("data/Atributos _org2011.csv", header=TRUE)

# Verificar dirección de las redes
is_directed <- list(
  Cooperacion = !isSymmetric(coop_matrix),
  Confianza    = !isSymmetric(trust_matrix),
  Recursos     = !isSymmetric(resources_matrix),
  Valores      = !isSymmetric(values_matrix),
  Parentesco   = !isSymmetric(parentesco_matrix)
)
cat("Verificación de direccionalidad de las redes:\n")
print(unlist(is_directed))

# Lista de matrices para análisis
matrices_list <- list(
  Cooperacion = coop_matrix,
  Confianza = trust_matrix,
  Recursos = resources_matrix,
  Valores = values_matrix,
  Parentesco = parentesco_matrix
)

# ---------------------------
# Estadísticas descriptivas de redes
# ---------------------------
cat("\nCalculando estadísticas descriptivas de las redes...\n")
descript_stats <- map_dfr(matrices_list, network_stats, .id = "Red")
descript_stats %>% 
  kbl(format = "latex", booktabs = TRUE, 
      caption = "Estadísticas descriptivas de las redes") %>% 
  kable_styling(latex_options = c("hold_position", "scale_down")) %>% 
  print()
# ---------------------------
# Análisis de Actores
# ---------------------------
cat("\nAnálisis de Actores...\n")
# Crear objeto network para la red de cooperación y agregar atributos
net_coop <- network(coop_matrix, directed = TRUE)
net_coop %v% "tipo"       <- attributes$Tipo_de_organización
net_coop %v% "ubicacion"  <- attributes$Ubicación
net_coop %v% "orientacion"<- attributes$Orientación

# Calcular medidas de centralidad
cent_df <- data.frame(
  Organizacion = attributes$Nombre,
  Tipo = attributes$Tipo_de_organización,
  Orientacion = attributes$Orientación,
  Degree = sna::degree(net_coop, gmode="digraph"),
  Indegree = sna::degree(net_coop, gmode="digraph", cmode="indegree"),
  Outdegree = sna::degree(net_coop, gmode="digraph", cmode="outdegree"),
  Betweenness = sna::betweenness(net_coop),
  Eigenvector = sna::evcent(net_coop)
)

# Top 10 actores según Indegree
cat("\nTop 10 actores según Indegree:\n")
cent_df %>% 
  arrange(desc(Indegree)) %>% 
  select(Organizacion, Tipo, Indegree) %>% 
  head(10) %>% 
  kbl(format = "latex", booktabs = TRUE, 
      caption = "") %>% 
  kable_styling(latex_options = c("hold_position", "scale_down")) %>% 
  print()
  

# Top 10 actores según Betweenness
cat("\nTop 10 actores según Betweenness:\n")
cent_df %>% 
  arrange(desc(Betweenness)) %>% 
  select(Organizacion, Tipo, Betweenness) %>% 
  head(10) %>% 
  kbl(format = "latex", booktabs = TRUE, 
      caption = "") %>% 
  kable_styling(latex_options = c("hold_position", "scale_down")) %>% 
  print()

# Centralidad promedio por tipo de organización
cat("\nCentralidad promedio por tipo de organización:\n")
cent_by_type <- cent_df %>% 
  group_by(Tipo) %>% 
  summarise(
    Mean_Indegree = mean(Indegree),
    Mean_Betweenness = mean(Betweenness),
    N = n()
  )
cent_by_type %>% 
  kbl(format = "latex", booktabs = TRUE, 
      caption = "") %>% 
  kable_styling(latex_options = c("hold_position", "scale_down")) %>% 
  print()

# ---------------------------
# Equivalencia estructural (CONCOR)
# ---------------------------
cat("\nAnálisis de Equivalencia Estructural (CONCOR)...\n")
# Se utiliza la red de cooperación para el análisis CONCOR
bloques <- concor(list(coop_matrix), cutoff = 0.8, nsplit = 2)
# Se espera que se identifiquen 4 bloques

for(i in 1:4) {
  cat(paste("\n\nBLOQUE", i, "\n"))
  cat(strrep("-", 15), "\n")
  organizaciones <- attributes[bloques == i, ]
  cat("Composición por Tipo:\n")
  print(table(organizaciones$Tipo_de_organización))
  cat("\nOrganizaciones específicas:\n")
  print(organizaciones$Nombre)
}

# Calcular y mostrar densidades entre bloques
block_density <- matrix(0, nrow = 4, ncol = 4)
for(i in 1:4){
  for(j in 1:4){
    nodes_i <- which(bloques == i)
    nodes_j <- which(bloques == j)
    if(length(nodes_i) > 0 & length(nodes_j) > 0){
      block_density[i,j] <- sum(coop_matrix[nodes_i, nodes_j]) / (length(nodes_i) * length(nodes_j))
    }
  }
}
cat("\n\nDensidad entre bloques:\n")
print(round(block_density, 3))

# Calcular cohesión de bloques y mostrar resultados
for(i in 1:4){
  cat(paste("\n\nCohesión del Bloque", i, "\n"))
  internal <- block_density[i,i]
  external <- mean(block_density[i, -i])
  ei_index <- (external - internal) / (external + internal)
  cat(paste("Densidad interna:", round(internal, 3), "\n"))
  cat(paste("Densidad externa media:", round(external, 3), "\n"))
  cat(paste("Índice E-I:", round(ei_index, 3), "\n"))
}

# Visualización: Heatmap de densidades y composición de bloques
cat("\nGenerando visualizaciones de bloques...\n")
# 1. DataFrame para heatmap de densidades entre bloques
density_df <- expand.grid(
  from = factor(1:4, labels = paste("Bloque", 1:4)),
  to = factor(1:4, labels = paste("Bloque", 1:4))
) %>% 
  mutate(density = as.vector(block_density))

# 2. DataFrame para la composición de bloques (se asume la composición definida)
composicion_bloques <- data.frame(
  Bloque = rep(1:4, c(18, 13, 14, 15)),
  Tipo = c(
    rep(c("Comite Vivienda", "Org. Cultural", "Club Deportivo", "Gob. Municipal"), c(9,4,4,1)), # Bloque 1
    rep(c("Org. Política Base", "Org. Cultural", "Otros"), c(1,8,4)),                             # Bloque 2
    rep(c("Org. Vecinal", "Org. Política Pobladores", "Org. Política Base"), c(3,6,5)),             # Bloque 3
    rep("Org. Vecinal", 15)                                                                         # Bloque 4
  )
)

# 3. DataFrame para métricas de cohesión
cohesion_df <- data.frame(
  Bloque = factor(1:4, labels = paste("Bloque", 1:4)),
  Densidad_interna = diag(block_density),
  Densidad_externa = sapply(1:4, function(i) mean(block_density[i, -i])),
  EI_index = sapply(1:4, function(i) {
    internal <- block_density[i,i]
    external <- mean(block_density[i, -i])
    (external - internal) / (external + internal)
  })
) %>% 
  pivot_longer(cols = -Bloque, names_to = "Metrica", values_to = "Valor")


# Custom theme for larger text, in English
my_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
}

# 1. Heatmap of block densities
p1 <- ggplot(density_df, aes(x = from, y = to, fill = density)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "viridis") +
  labs(
    title = "Density Between Blocks",
    x = "Origin Block",
    y = "Destination Block",
    fill = "Density"
  ) +
  my_theme() +
  theme(axis.text.x = element_text( vjust = 1, hjust = 1))

# 2. Block composition
p2 <- ggplot(composicion_bloques, aes(x = factor(Bloque), fill = Tipo)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Block Composition",
    x = "Block",
    y = "Proportion",
    fill = "Organization Type"
  ) +
  my_theme()

# 3. Cohesion metrics per block
p3 <- ggplot(cohesion_df, aes(x = Bloque, y = Valor, fill = Metrica)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Cohesion Metrics by Block",
    x = "Block",
    y = "Value",
    fill = "Metric"
  ) +
  my_theme()

# Convert each ggplot to plotly for interactivity
p1_plotly <- ggplotly(p1)
p2_plotly <- ggplotly(p2)
p3_plotly <- ggplotly(p3)

# Print or display the interactive plots
print(p1_plotly)
print(p2_plotly)
print(p3_plotly)

# ---------------------------
# Análisis de Red Multiplex
# ---------------------------
cat("\nConstruyendo red multiplex...\n")
multiplex_matrix <- coop_matrix * 0.3 + 
  trust_matrix * 0.3 + 
  resources_matrix * 0.2 + 
  values_matrix * 0.2

g <- graph_from_adjacency_matrix(multiplex_matrix, mode = "directed", weighted = TRUE)
V(g)$degree      <- igraph::degree(g)
V(g)$betweenness <- igraph::betweenness(g)
V(g)$eigenvector <- eigen_centrality(g)$vector
V(g)$block       <- attributes$X1PosciónCONCOR  # Asegúrate que este campo exista en tus atributos
V(g)$type        <- attributes$Tipo_de_organización
V(g)$orientation <- attributes$Orientación
V(g)$nombre      <- attributes$Nombre

# Agregar atributo para el tipo de vínculo según pertenencia a bloques
edge_list <- as_edgelist(g)
E(g)$block_tie <- ifelse(
  V(g)$block[edge_list[,1]] == V(g)$block[edge_list[,2]], 
  paste0("Block", V(g)$block[edge_list[,1]]), 
  "between"
)

# Graficar la red multiplex con ggraph
p_net <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(alpha = ifelse(V(g)$block[from] == V(g)$block[to], 0.6, 0.2)),
                 color = "black",
                 arrow = arrow(length = unit(2, 'mm'))) +
  geom_node_point(aes(size = betweenness,
                      color = as.factor(block)),
                  alpha = 0.7) +
  geom_node_label(aes(label = ifelse(betweenness > quantile(betweenness, 0.92), nombre, "")),
                  repel = TRUE,
                  size = 3) +
  scale_color_viridis_d(option = "D",
                        name = "Bloque",
                        labels = c("Político", "JJ.VV", "Base", "Cultural")) +
  scale_size_continuous(name = "Intermediación",
                        range = c(2, 8),
                        breaks = c(0, 0.025, 0.05, 0.075, 0.1, 0.125),
                        labels = c("0.000", "0.025", "0.050", "0.075", "0.100", "0.125")) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         size = guide_legend(override.aes = list(alpha = 1)),
         edge_alpha = "none") +
  labs(title = "Multiplex relations network") +
  theme_graph() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

print(p_net)

cat("\nAnálisis completado.\n")



