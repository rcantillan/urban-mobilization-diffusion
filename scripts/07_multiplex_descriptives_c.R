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
# 1. Funciones de utilidad
# ============================================================================

#' Calcula métricas básicas para una capa de red (usa formato igraph)
#' @param g Objeto igraph
#' @return Lista con métricas estructurales básicas
calculate_layer_metrics <- function(g) {
  list(
    density       = edge_density(g),
    avg_degree    = mean(degree(g)),
    reciprocity   = reciprocity(g),
    clustering    = transitivity(g),
    components    = components(g)$no,
    mean_distance = mean_distance(g, directed = FALSE)
  )
}

#' Analiza relaciones entre pares de capas (usa matrices dispersas para calcular solapamientos)
#' @param mat1,mat2 Matrices de adyacencia (sparse matrix recomendadas)
calculate_layer_relations <- function(mat1, mat2) {
  # Convertir a lógicas para encontrar intersecciones y uniones
  logical1 <- mat1 > 0
  logical2 <- mat2 > 0
  
  overlap       <- sum(logical1 & logical2)
  union_size    <- sum(logical1 | logical2)
  jaccard       <- ifelse(union_size > 0, overlap / union_size, 0)
  
  # Correlación estructural entre capas
  v1          <- as.numeric(mat1)
  v2          <- as.numeric(mat2)
  correlation <- cor(v1, v2)
  
  list(
    overlap_count = overlap,
    jaccard       = jaccard,
    correlation   = correlation
  )
}

#' Calcula centralidad multiplexada
#' @param mats Lista de matrices de adyacencia (en formato igraph o matrix)
#' @return Data frame con medidas de centralidad por actor
calculate_multiplex_centrality <- function(mats) {
  n_actors <- nrow(mats[[1]])
  
  degree_mat   <- matrix(0, nrow = n_actors, ncol = length(mats))
  between_mat  <- matrix(0, nrow = n_actors, ncol = length(mats))
  
  for (i in seq_along(mats)) {
    g <- if (inherits(mats[[i]], "igraph")) {
      mats[[i]]
    } else {
      graph_from_adjacency_matrix(mats[[i]])
    }
    degree_mat[, i]  <- degree(g)
    between_mat[, i] <- betweenness(g)
    
    rm(g)
    gc()
  }
  
  data.frame(
    actor        = seq_len(n_actors),
    degree_mean  = rowMeans(degree_mat),
    degree_sd    = apply(degree_mat, 1, sd),
    between_mean = rowMeans(between_mat),
    between_sd   = apply(between_mat, 1, sd),
    n_layers     = rowSums(degree_mat > 0)
  )
}

# ============================================================================
# Para TABLA COMPARATIVA de comunidades (al estilo del manual)
# ============================================================================

#' Calcula estadísticas de comunidad inspiradas en el manual "Analysis of Multiplex Social Networks with R"
#' @param comm_result Data frame con columnas: actor, layer, cid (comunidad)
#' @param net Objeto multinet
#' @return Data frame con num_communities, avg_size, prop_vertices, prop_actors, actor_overl
community_stats <- function(comm_result, net) {
  # Número de comunidades
  num_com <- length(unique(comm_result$cid))
  
  # Tamaño de cada comunidad (en términos de vértices actor-layer)
  sizes <- table(comm_result$cid)
  avg_size <- mean(sizes)
  
  # Número total de vértices en la red (actor-capa)
  total_vertices <- nrow(vertices_ml(net))
  covered_verts  <- nrow(comm_result)  # filas = vértices en alguna comunidad
  prop_vertices  <- round(covered_verts / total_vertices, 3)
  
  # Proporción de actores cubiertos
  all_actors      <- unique(actors_ml(net))
  assigned_actors <- unique(comm_result$actor)
  prop_actors     <- round(length(intersect(all_actors, assigned_actors)) / length(all_actors), 3)
  
  # Grado de solapamiento a nivel de actor
  # (# pares actor-comunidad) / (# actores asignados)
  ac_pairs <- nrow(distinct(comm_result, actor, cid))
  actor_overl <- if (length(assigned_actors) > 0) {
    round(ac_pairs / length(assigned_actors), 3)
  } else {
    0
  }
  
  data.frame(
    num_communities = num_com,
    avg_size        = round(avg_size, 2),
    prop_vertices   = prop_vertices,
    prop_actors     = prop_actors,
    actor_overl     = actor_overl
  )
}

# ============================================================================
# 2. Lectura y preparación de datos (matrices dispersas)
# ============================================================================

mats_paths <- list(
  Cooperation = "data/CoopNet.csv",
  Trust       = "data/ConfianzaNet.csv",
  Resources   = "data/RecursosNet.csv", 
  Values      = "data/ValoresNet.csv",
  Kinship     = "data/ParentescoNet.csv"
)

# Función para leer matrices de adyacencia en formato disperso
read_network_matrix_sparse <- function(path) {
  if (!file.exists(path)) stop(paste("Archivo no encontrado:", path))
  mat <- as.matrix(read.csv(path, header = FALSE))
  return(Matrix(mat, sparse = TRUE))
}

cat("Leyendo y convirtiendo matrices de adyacencia a formato disperso...\n")
m_layers <- lapply(mats_paths, read_network_matrix_sparse)
gc()

# ============================================================================
# 3. Análisis por capa
# ============================================================================

layer_metrics_list <- list()
for (layer_name in names(m_layers)) {
  g_temp <- graph_from_adjacency_matrix(as.matrix(m_layers[[layer_name]]))
  layer_metrics_list[[layer_name]] <- calculate_layer_metrics(g_temp)
  rm(g_temp)
  gc()
}

layer_metrics_df <- do.call(rbind, lapply(names(layer_metrics_list), function(layer) {
  metrics <- layer_metrics_list[[layer]]
  data.frame(
    Layer         = layer,
    Density       = round(metrics$density, 3),
    AvgDegree     = round(metrics$avg_degree, 2),
    Reciprocity   = round(metrics$reciprocity, 3),
    Clustering    = round(metrics$clustering, 3),
    Components    = metrics$components,
    MeanDistance  = round(metrics$mean_distance, 2)
  )
}))

# Mostrar en LaTeX (o consola)
kbl(layer_metrics_df, 
    format   = "latex", 
    caption  = "Métricas estructurales por capa de red",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 4. Análisis de interrelaciones entre capas
# ============================================================================

layer_pairs <- combn(names(m_layers), 2, simplify = FALSE)
layer_relations <- list()

for (pair in layer_pairs) {
  rel <- calculate_layer_relations(m_layers[[pair[1]]], m_layers[[pair[2]]])
  layer_relations[[paste(pair, collapse = "_")]] <- data.frame(
    Layer1      = pair[1],
    Layer2      = pair[2],
    Overlap     = rel$overlap_count,
    Jaccard     = round(rel$jaccard, 3),
    Correlation = round(rel$correlation, 3)
  )
}

layer_relations_df <- do.call(rbind, layer_relations)

kbl(layer_relations_df,
    format   = "latex",
    caption  = "Relaciones estructurales entre pares de capas",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 5. Centralidad multiplexada
# ============================================================================

centrality_df <- calculate_multiplex_centrality(m_layers)

actor_types <- kmeans(
  scale(centrality_df[, c("degree_mean", "between_mean", "n_layers")]), 
  centers = 4
)
centrality_df$type <- factor(actor_types$cluster)

actor_summary <- centrality_df %>%
  group_by(type) %>%
  summarise(
    n            = n(),
    degree_mean  = mean(degree_mean),
    between_mean = mean(between_mean),
    avg_layers   = mean(n_layers)
  )

kbl(actor_summary,
    format   = "latex",
    caption  = "Tipos de actores según centralidad multiplexada",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 6. Guardar resultados en archivos
# ============================================================================

if (!dir.exists("output")) dir.create("output")
write.csv(layer_metrics_df,     "output/layer_metrics.csv",    row.names = FALSE)
write.csv(layer_relations_df,   "output/layer_relations.csv",  row.names = FALSE)
write.csv(centrality_df,        "output/actor_centrality.csv", row.names = FALSE)

# ============================================================================
# 7. Construcción del objeto multinet
# ============================================================================
gc()
cat("Creando objeto multinet...\n")
net <- ml_empty()
add_layers_ml(net, names(m_layers))

cat("Agregando vértices y aristas a cada capa...\n")
for (layer_name in names(m_layers)) {
  mat_layer <- m_layers[[layer_name]]
  
  # Convertir la matriz dispersa a indices (row, col) donde hay enlaces
  edges <- summary(mat_layer)  # indices de celdas != 0 para dgCMatrix
  
  if (nrow(edges) > 0) {
    edges_df <- data.frame(
      actors_from = edges$i,
      layers_from = layer_name,
      actors_to   = edges$j,
      layers_to   = layer_name
    )
    add_edges_ml(net, edges_df)
  }
  rm(edges, mat_layer)
  gc()
}

# ============================================================================
# 8. Detección de comunidades con ABACUS, Clique, Louvain e Infomap
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
# 1. Funciones de utilidad
# ============================================================================

# Métricas por capa (formato igraph)
calculate_layer_metrics <- function(g) {
  list(
    density       = edge_density(g),
    avg_degree    = mean(degree(g)),
    reciprocity   = reciprocity(g),
    clustering    = transitivity(g),
    components    = components(g)$no,
    mean_distance = mean_distance(g, directed = FALSE)
  )
}

# Relaciones entre pares de capas (usa sparse matrix)
calculate_layer_relations <- function(mat1, mat2) {
  logical1 <- mat1 > 0
  logical2 <- mat2 > 0
  
  overlap    <- sum(logical1 & logical2)
  union_size <- sum(logical1 | logical2)
  jaccard    <- ifelse(union_size > 0, overlap / union_size, 0)
  
  v1          <- as.numeric(mat1)
  v2          <- as.numeric(mat2)
  correlation <- cor(v1, v2)
  
  list(
    overlap_count = overlap,
    jaccard       = jaccard,
    correlation   = correlation
  )
}

# Centralidad "multiplex"
calculate_multiplex_centrality <- function(mats) {
  n_actors <- nrow(mats[[1]])
  
  degree_mat   <- matrix(0, nrow = n_actors, ncol = length(mats))
  between_mat  <- matrix(0, nrow = n_actors, ncol = length(mats))
  
  for (i in seq_along(mats)) {
    g <- if (inherits(mats[[i]], "igraph")) {
      mats[[i]]
    } else {
      graph_from_adjacency_matrix(mats[[i]])
    }
    degree_mat[, i]  <- degree(g)
    between_mat[, i] <- betweenness(g)
    
    rm(g)
    gc()
  }
  
  data.frame(
    actor        = seq_len(n_actors),
    degree_mean  = rowMeans(degree_mat),
    degree_sd    = apply(degree_mat, 1, sd),
    between_mean = rowMeans(between_mat),
    between_sd   = apply(between_mat, 1, sd),
    n_layers     = rowSums(degree_mat > 0)
  )
}

# Estadísticas de comunidad (inspirado en el manual "Analysis of Multiplex Social Networks with R")
community_stats <- function(comm_result, net) {
  # Si está vacío o mal formado, devuelves algo minimal
  if (is.null(comm_result) || nrow(comm_result) == 0 || !"cid" %in% names(comm_result)) {
    return(data.frame(
      num_communities = 0,
      avg_size        = NA,
      prop_vertices   = 0,
      prop_actors     = 0,
      actor_overl     = 0
    ))
  }
  
  num_com <- length(unique(comm_result$cid))
  sizes   <- table(comm_result$cid)
  avg_size <- mean(sizes)
  
  total_vertices <- nrow(vertices_ml(net))
  covered_verts  <- nrow(comm_result)
  prop_vertices  <- round(covered_verts / total_vertices, 3)
  
  all_actors      <- unique(actors_ml(net))
  assigned_actors <- unique(comm_result$actor)
  prop_actors     <- round(length(intersect(all_actors, assigned_actors)) / length(all_actors), 3)
  
  ac_pairs     <- nrow(dplyr::distinct(comm_result, actor, cid))
  actor_overl  <- if (length(assigned_actors) > 0) {
    round(ac_pairs / length(assigned_actors), 3)
  } else {
    0
  }
  
  data.frame(
    num_communities = num_com,
    avg_size        = round(avg_size, 2),
    prop_vertices   = prop_vertices,
    prop_actors     = prop_actors,
    actor_overl     = actor_overl
  )
}

# ============================================================================
# 2. Lectura y preparación de datos (matrices dispersas)
# ============================================================================

mats_paths <- list(
  Cooperation = "data/CoopNet.csv",
  Trust       = "data/ConfianzaNet.csv",
  Resources   = "data/RecursosNet.csv", 
  Values      = "data/ValoresNet.csv",
  Kinship     = "data/ParentescoNet.csv"
)

read_network_matrix_sparse <- function(path) {
  if (!file.exists(path)) stop(paste("Archivo no encontrado:", path))
  mat <- as.matrix(read.csv(path, header = FALSE))
  return(Matrix(mat, sparse = TRUE))
}

cat("Leyendo y convirtiendo matrices de adyacencia a formato disperso...\n")
m_layers <- lapply(mats_paths, read_network_matrix_sparse)
gc()

# ============================================================================
# 3. Análisis por capa
# ============================================================================

layer_metrics_list <- list()
for (layer_name in names(m_layers)) {
  g_temp <- graph_from_adjacency_matrix(as.matrix(m_layers[[layer_name]]))
  layer_metrics_list[[layer_name]] <- calculate_layer_metrics(g_temp)
  rm(g_temp)
  gc()
}

layer_metrics_df <- do.call(rbind, lapply(names(layer_metrics_list), function(layer) {
  metrics <- layer_metrics_list[[layer]]
  data.frame(
    Layer         = layer,
    Density       = round(metrics$density, 3),
    AvgDegree     = round(metrics$avg_degree, 2),
    Reciprocity   = round(metrics$reciprocity, 3),
    Clustering    = round(metrics$clustering, 3),
    Components    = metrics$components,
    MeanDistance  = round(metrics$mean_distance, 2)
  )
}))

kbl(layer_metrics_df, 
    format   = "latex", 
    caption  = "Métricas estructurales por capa de red",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 4. Análisis de interrelaciones entre capas
# ============================================================================

layer_pairs <- combn(names(m_layers), 2, simplify = FALSE)
layer_relations <- list()

for (pair in layer_pairs) {
  rel <- calculate_layer_relations(m_layers[[pair[1]]], m_layers[[pair[2]]])
  layer_relations[[paste(pair, collapse = "_")]] <- data.frame(
    Layer1      = pair[1],
    Layer2      = pair[2],
    Overlap     = rel$overlap_count,
    Jaccard     = round(rel$jaccard, 3),
    Correlation = round(rel$correlation, 3)
  )
}

layer_relations_df <- do.call(rbind, layer_relations)

kbl(layer_relations_df,
    format   = "latex",
    caption  = "Relaciones estructurales entre pares de capas",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 5. Centralidad multiplexada
# ============================================================================

centrality_df <- calculate_multiplex_centrality(m_layers)

actor_types <- kmeans(
  scale(centrality_df[, c("degree_mean", "between_mean", "n_layers")]), 
  centers = 4
)
centrality_df$type <- factor(actor_types$cluster)

actor_summary <- centrality_df %>%
  group_by(type) %>%
  summarise(
    n            = n(),
    degree_mean  = mean(degree_mean),
    between_mean = mean(between_mean),
    avg_layers   = mean(n_layers)
  )

kbl(actor_summary,
    format   = "latex",
    caption  = "Tipos de actores según centralidad multiplexada",
    booktabs = TRUE) %>%
  kable_styling(latex_options = "HOLD_position")

# ============================================================================
# 6. Guardar resultados en archivos
# ============================================================================

if (!dir.exists("output")) dir.create("output")
write.csv(layer_metrics_df,     "output/layer_metrics.csv",    row.names = FALSE)
write.csv(layer_relations_df,   "output/layer_relations.csv",  row.names = FALSE)
write.csv(centrality_df,        "output/actor_centrality.csv", row.names = FALSE)

# ============================================================================
# 7. Construcción del objeto multinet
# ============================================================================
