# ============================================================================
# Análisis de Redes Múltiples (Multiplex)
# Basado en: Magnani & Rossi (2011) y otros trabajos académicos
# ============================================================================
gc()
suppressMessages({
  library(tidyverse)
  library(igraph)
  library(Matrix)
  library(knitr)
  library(kableExtra)
  library(multinet)
})

# ============================================================================
# 1. Funciones de utilidad para métricas multiplexadas
# ============================================================================

#' Calcula métricas básicas para una capa de red
#' @param mat Matriz de adyacencia
#' @return Lista con métricas estructurales básicas
calculate_layer_metrics <- function(mat) {
  g <- graph_from_adjacency_matrix(mat)
  list(
    density = edge_density(g),
    avg_degree = mean(degree(g)),
    reciprocity = reciprocity(g),
    clustering = transitivity(g),
    components = components(g)$no,
    mean_distance = mean_distance(g, directed=FALSE)
  )
}

#' Analiza relaciones entre pares de capas
#' @param mat1,mat2 Par de matrices de adyacencia a comparar
calculate_layer_relations <- function(mat1, mat2) {
  # Calculamos solapamiento usando índice de Jaccard
  overlap <- sum((mat1 > 0) & (mat2 > 0))
  union_size <- sum((mat1 > 0) | (mat2 > 0))
  jaccard <- ifelse(union_size > 0, overlap/union_size, 0)
  
  # Correlación estructural entre capas
  correlation <- cor(as.vector(mat1), as.vector(mat2))
  
  list(
    overlap_count = overlap,
    jaccard = jaccard, 
    correlation = correlation
  )
}

#' Calcula centralidad multiplexada
#' @param mats Lista de matrices de adyacencia
#' @return Data frame con medidas de centralidad por actor
calculate_multiplex_centrality <- function(mats) {
  n_actors <- nrow(mats[[1]])
  
  # Inicializamos matrices para almacenar resultados
  degree_mat <- matrix(0, nrow=n_actors, ncol=length(mats))
  between_mat <- matrix(0, nrow=n_actors, ncol=length(mats))
  
  # Calculamos centralidad por capa
  for(i in seq_along(mats)) {
    g <- graph_from_adjacency_matrix(mats[[i]])
    degree_mat[,i] <- degree(g)
    between_mat[,i] <- betweenness(g)
  }
  
  # Calculamos medidas agregadas
  data.frame(
    actor = 1:n_actors,
    degree_mean = rowMeans(degree_mat),
    degree_sd = apply(degree_mat, 1, sd),
    between_mean = rowMeans(between_mat),
    between_sd = apply(between_mat, 1, sd),
    n_layers = rowSums(degree_mat > 0)
  )
}

# ============================================================================
# 2. Lectura y preparación de datos
# ============================================================================

# Definimos rutas a los archivos de red
mats_paths <- list(
  Cooperation = "data/CoopNet.csv",
  Trust = "data/ConfianzaNet.csv",
  Resources = "data/RecursosNet.csv", 
  Values = "data/ValoresNet.csv",
  Kinship = "data/ParentescoNet.csv"
)

# Función para leer matrices de adyacencia
read_network_matrix <- function(path) {
  if(!file.exists(path)) stop(paste("Archivo no encontrado:", path))
  as.matrix(read.csv(path, header=FALSE))
}

cat("Leyendo matrices de adyacencia...\n")
m_layers <- lapply(mats_paths, read_network_matrix)

# ============================================================================
# 3. Análisis de estructura por capa
# ============================================================================

# Calculamos métricas básicas por capa
layer_metrics <- lapply(m_layers, calculate_layer_metrics)

# Creamos tabla de métricas por capa
layer_metrics_df <- do.call(rbind, lapply(names(layer_metrics), function(layer) {
  metrics <- layer_metrics[[layer]]
  data.frame(
    Layer = layer,
    Density = round(metrics$density, 3),
    AvgDegree = round(metrics$avg_degree, 2),
    Reciprocity = round(metrics$reciprocity, 3),
    Clustering = round(metrics$clustering, 3),
    Components = metrics$components,
    MeanDistance = round(metrics$mean_distance, 2)
  )
}))

# Imprimimos tabla en LaTeX
kbl(layer_metrics_df, 
    format = "latex", 
    caption = "Métricas estructurales por capa de red",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 4. Análisis de interrelaciones entre capas
# ============================================================================

# Calculamos relaciones entre todos los pares de capas
layer_pairs <- combn(names(m_layers), 2, simplify = FALSE)
layer_relations <- lapply(layer_pairs, function(pair) {
  rel <- calculate_layer_relations(m_layers[[pair[1]]], m_layers[[pair[2]]])
  data.frame(
    Layer1 = pair[1],
    Layer2 = pair[2],
    Overlap = rel$overlap_count,
    Jaccard = round(rel$jaccard, 3),
    Correlation = round(rel$correlation, 3)
  )
})

# Creamos tabla de relaciones entre capas
layer_relations_df <- do.call(rbind, layer_relations)

kbl(layer_relations_df,
    format = "latex",
    caption = "Relaciones estructurales entre pares de capas",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 5. Análisis de centralidad multiplexada
# ============================================================================

# Calculamos medidas de centralidad multiplexada
centrality_df <- calculate_multiplex_centrality(m_layers)

# Identificamos tipos de actores según patrones de centralidad
actor_types <- kmeans(
  scale(centrality_df[,c("degree_mean", "between_mean", "n_layers")]), 
  centers = 4
)
centrality_df$type <- factor(actor_types$cluster)

# Creamos tabla resumen por tipo de actor
actor_summary <- centrality_df %>%
  group_by(type) %>%
  summarise(
    n = n(),
    degree_mean = mean(degree_mean),
    between_mean = mean(between_mean),
    avg_layers = mean(n_layers)
  )

kbl(actor_summary,
    format = "latex",
    caption = "Tipos de actores según centralidad multiplexada",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 6. Guardado de resultados
# ============================================================================

# Guardamos tablas principales en archivos
write.csv(layer_metrics_df, "output/layer_metrics.csv", row.names=FALSE)
write.csv(layer_relations_df, "output/layer_relations.csv", row.names=FALSE)
write.csv(centrality_df, "output/actor_centrality.csv", row.names=FALSE)

# ============================================================================
# Análisis de Comunidades en Redes Multiplex
# Basado en el manual oficial de multinet v4.2.1
# ============================================================================

#' Análisis completo de comunidades usando todos los métodos de multinet
#' @param net Red multiplex
#' @return Lista con resultados de diferentes algoritmos
analyze_multiplex_communities <- function(net) {
  results <- list()
  
  # 1. Método Louvain Generalizado
  # Este método optimiza la modularidad considerando tanto conexiones dentro de capas
  # como la consistencia de comunidades entre capas
  cat("Analizando comunidades con método Louvain...\n")
  tryCatch({
    results$louvain <- glouvain_ml(net)
    cat("Número de comunidades (Louvain):", length(unique(results$louvain$cid)), "\n")
  }, error = function(e) {
    cat("Error en método Louvain:", e$message, "\n")
  })
  
  # 2. Método ABACUS
  # Para comunidades superpuestas usando minería de patrones frecuentes
  #cat("\nAnalizando comunidades con ABACUS...\n")
  #tryCatch({
  #  results$abacus <- abacus_ml(net, 4, 2)  # k=4 nodos, m=2 capas mínimas
  #  cat("Número de comunidades (ABACUS):", length(unique(results$abacus$cid)), "\n")
  #}, error = function(e) {
  #  cat("Error en método ABACUS:", e$message, "\n")
  #})
  
  # 3. Método Clique Percolation
  # Para detectar comunidades superpuestas basadas en cliques
  cat("\nAnalizando comunidades con Clique Percolation...\n")
  tryCatch({
    results$clique <- clique_percolation_ml(net, 4, 2)  # k=4 tamaño clique, m=2 capas
    cat("Número de comunidades (Clique):", length(unique(results$clique$cid)), "\n")
  }, error = function(e) {
    cat("Error en método Clique Percolation:", e$message, "\n")
  })
  
  # 4. Método Infomap
  # Basado en la compresión de la descripción de caminatas aleatorias
  cat("\nAnalizando comunidades con Infomap...\n")
  tryCatch({
    results$infomap <- infomap_ml(net)
    cat("Número de comunidades (Infomap):", length(unique(results$infomap$cid)), "\n")
  }, error = function(e) {
    cat("Error en método Infomap:", e$message, "\n")
  })
  
  # Comparar resultados de los métodos exitosos
  successful_methods <- !sapply(results, is.null)
  if(sum(successful_methods) > 0) {
    stats_df <- data.frame(
      method = names(results)[successful_methods],
      n_communities = sapply(results[successful_methods], function(x) {
        length(unique(x$cid))
      }),
      avg_size = sapply(results[successful_methods], function(x) {
        mean(table(x$cid))
      }),
      clust_vertices = sapply(results[successful_methods], function(x) {
        length(unique(x$actor))
      }),
      clust_actors = sapply(results[successful_methods], function(x) {
        n_distinct(x$actor[!duplicated(paste(x$actor, x$cid))])
      })
    )
    results$comparison <- stats_df
  }
  
  return(results)
}

#' Visualiza las comunidades detectadas usando el layout multiforce
#' @param net Red multiplex
#' @param communities Resultado de detección de comunidades
plot_communities_ml <- function(net, communities, method_name) {
  # Calcular layout considerando la estructura multicapa
  l <- layout_multiforce_ml(net, w_inter = 0.1, gravity = 1)
  
  # Generar paleta de colores para las comunidades
  n_comm <- length(unique(communities$cid))
  comm_colors <- rainbow(n_comm)
  
  # Visualizar con configuración del manual
  plot(net, 
       layout = l,
       grid = c(2, 3),
       vertex.color = comm_colors[as.factor(communities$cid)],
       vertex.labels = "",
       vertex.size = 3,
       main = paste("Comunidades -", method_name))
}

# Ejecutar análisis
cat("Cargando datos...\n")
mats <- lapply(mats_paths, read_network_matrix)

# Crear objeto multinet
net <- ml_empty()
add_layers_ml(net, names(mats))

# Agregar vértices y aristas por capa
for(layer_name in names(mats)) {
  cat("Procesando capa", layer_name, "...\n")
  mat <- mats[[layer_name]]
  edges <- which(mat > 0, arr.ind = TRUE)
  
  if(nrow(edges) > 0) {
    edges_df <- data.frame(
      actors_from = edges[,1],
      layers_from = layer_name,
      actors_to = edges[,2],
      layers_to = layer_name
    )
    add_edges_ml(net, edges_df)
  }
  gc()  # Limpiar memoria después de cada capa
}

# Ejecutar análisis de comunidades
results <- analyze_multiplex_communities(net)

# Mostrar tabla comparativa
if(!is.null(results$comparison)) {
  cat("\nComparación de métodos de detección de comunidades:\n")
  print(results$comparison)
  
  # Visualizar comunidades para cada método exitoso
  for(method in names(results)[!names(results) %in% "comparison"]) {
    if(!is.null(results[[method]])) {
      cat("\nVisualizando comunidades detectadas con", method, "...\n")
      plot_communities_ml(net, results[[method]], method)
    }
  }
} else {
  cat("\nNo se pudieron completar los análisis de comunidades\n")
}

# Análisis de robustez variando omega (solo para Louvain)
if(!is.null(results$louvain)) {
  cat("\nAnalizando sensibilidad al parámetro omega...\n")
  omega_values <- seq(0, 1, by = 0.2)
  omega_results <- lapply(omega_values, function(w) {
    comm <- glouvain_ml(net, omega = w)
    data.frame(
      omega = w,
      n_communities = length(unique(comm$cid)),
      avg_size = mean(table(comm$cid))
    )
  })
  omega_comparison <- do.call(rbind, omega_results)
  
  # Visualizar efecto de omega
  print(ggplot(omega_comparison, aes(x = omega)) +
          geom_line(aes(y = n_communities, color = "Número de comunidades")) +
          geom_line(aes(y = avg_size, color = "Tamaño promedio")) +
          labs(title = "Efecto del parámetro omega en la detección de comunidades",
               x = "Omega",
               y = "Valor",
               color = "Métrica") +
          theme_minimal())
}