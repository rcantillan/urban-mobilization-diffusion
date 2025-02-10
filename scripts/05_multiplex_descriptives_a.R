suppressMessages({
  library(tidyverse)
  library(igraph)
  library(Matrix)
  library(knitr)
  library(kableExtra)
  library(multinet)  # Asegúrate de tener la versión >= 4.0
})

###############################################################################
# (1) Lectura de datos desde csv, creando las matrices de adyacencia
###############################################################################
mats_paths <- list(
  Cooperation = "data/CoopNet.csv",
  Trust       = "data/ConfianzaNet.csv",
  Resources   = "data/RecursosNet.csv",
  Values      = "data/ValoresNet.csv",
  Kinship     = "data/ParentescoNet.csv"
)

read_network_matrix <- function(path) {
  if(!file.exists(path)) stop(paste("File not found:", path))
  as.matrix(read.csv(path, header=FALSE))
}

cat("Leyendo las matrices de adyacencia...\n")
m_layers <- lapply(mats_paths, read_network_matrix)

###############################################################################
# (2) Remover nodos inactivos en TODAS las capas
###############################################################################
cat("\nRemoviendo nodos con grado global cero...\n")

n_nodes <- nrow(m_layers[[1]])
degrees_total <- rep(0, n_nodes)

for(k in seq_along(m_layers)) {
  g_k <- graph_from_adjacency_matrix(m_layers[[k]], mode="directed", diag=FALSE)
  deg_k <- degree(g_k, mode="total")
  degrees_total <- degrees_total + deg_k
}

inactive_idx <- which(degrees_total == 0)
cat("Número de nodos con grado cero (global):", length(inactive_idx), "\n")

if(length(inactive_idx) > 0){
  for(k in seq_along(m_layers)){
    mat_k <- m_layers[[k]]
    mat_k <- mat_k[-inactive_idx, -inactive_idx, drop=FALSE]
    m_layers[[k]] <- mat_k
  }
  cat("Nodos inactivos removidos en todas las matrices.\n")
} else {
  cat("No se removieron nodos; todos están activos al menos en una capa.\n")
}

n_nodes <- nrow(m_layers[[1]])  # nueva dimensión

###############################################################################
# (3) Estadísticas descriptivas (Overlap, Participation, summary multinet...)
###############################################################################
cat("\n*** Estadísticas descriptivas: Overlap y Participation ***\n")

## (3a) Crear un array 3D con las capas
arr_3d <- array(0, dim=c(n_nodes, n_nodes, length(m_layers)),
                dimnames=list(NULL, NULL, names(m_layers)))
i_layer <- 1
for(nm in names(m_layers)){
  arr_3d[,,i_layer] <- m_layers[[nm]]
  i_layer <- i_layer + 1
}

## Overlap: cuántas capas contiene la arista (i->j)
overlap_matrix <- matrix(0, nrow=n_nodes, ncol=n_nodes)
for(r in 1:n_nodes){
  for(c in 1:n_nodes){
    overlap_matrix[r,c] <- sum(arr_3d[r,c,] != 0)
  }
}
overlap_vals <- overlap_matrix[overlap_matrix>0]
tab_overlap  <- table(overlap_vals)
perc_overlap <- round(100*tab_overlap / sum(tab_overlap),2)

cat("\nDistribución de overlap (cuántas capas comparten la arista):\n")
print(tab_overlap)
cat("\nDistribución porcentual del overlap:\n")
print(perc_overlap)

## (3b) Participación nodal: cuántas capas está activo cada nodo
count_node_multiplex <- function(node_id) {
  out <- 0
  for(k in seq_along(m_layers)){
    mat_k <- m_layers[[k]]
    out_deg <- sum(mat_k[node_id, ] != 0)
    in_deg  <- sum(mat_k[, node_id] != 0)
    if((out_deg + in_deg) > 0) out <- out + 1
  }
  return(out)
}
part_vec <- sapply(1:n_nodes, count_node_multiplex)
cat("\nParticipación de nodos en capas:\n")
cat("Media:", mean(part_vec), "SD:", sd(part_vec),
    "Mín:", min(part_vec), "Máx:", max(part_vec), "\n")

df_participacion <- data.frame(
  Nodo=1:min(10,n_nodes),
  Participacion=part_vec[1:min(10,n_nodes)]
)

kbl(df_participacion, format="latex", booktabs=TRUE,
    caption="Ejemplo de participación nodal (primeros 10 nodos)") %>%
  kable_styling(latex_options="HOLD_position") %>%
  print()

## (3c) Crear objeto multinet (para funciones avanzadas)
cat("\nCreando objeto multinet para análisis avanzado...\n")
net_ml <- ml_empty("MultiplexData")
for(nm in names(m_layers)){
  g_k <- graph_from_adjacency_matrix(
    m_layers[[nm]],
    mode="directed",
    diag=FALSE
  )
  add_igraph_layer_ml(net_ml, g_k, nm)
}

cat("\nResumen básico del objeto multinet:\n")
print(net_ml)

## Summary integrado de multinet
cat("\n*** Summary de la red multiplex (por capa y aplanado) ***\n")
df_sum <- summary(net_ml)
print(df_sum)

kbl(df_sum, format="latex", booktabs=TRUE,
    caption="Resumen (summary) de la red multiplex por capa y aplanado (_flat_)") %>%
  kable_styling(latex_options="HOLD_position") %>%
  print()

## (3d) Comparación de capas con layer_comparison_ml
cat("\n*** Comparación de capas (layer_comparison_ml) ***\n")
methods_demo <- c("jaccard.actors", "pearson.degree", "jaccard.edges")
list_comparisons <- list()

for(method_ in methods_demo){
  tmp_result <- layer_comparison_ml(net_ml, method=method_)
  tmp_as_matrix <- try(as.matrix(tmp_result), silent=TRUE)
  
  list_comparisons[[method_]] <- tmp_result
  
  cat("\nMétodo:", method_, "\n")
  print(tmp_result)
  
  # Si la salida se comporta como matriz bidimensional => imprimimos tabla en LaTeX
  if(!inherits(tmp_as_matrix,"try-error") &&
     is.matrix(tmp_as_matrix) &&
     nrow(tmp_as_matrix)>1 && ncol(tmp_as_matrix)>1) {
    
    tmp_rounded <- round(tmp_as_matrix, 3)
    tmp_df <- as.data.frame(as.table(tmp_rounded))
    colnames(tmp_df) <- c("Capa1","Capa2","Valor")
    
    kbl(tmp_df, format="latex", booktabs=TRUE,
        caption=paste("Comparación de capas con método:", method_)) %>%
      kable_styling(latex_options="HOLD_position") %>%
      print()
  } else {
    cat("Aviso: la salida no es una matriz bidimensional, no se creará tabla.\n")
  }
}

###############################################################################
# (4a) Detección de comunidades (Agregación ponderada)
###############################################################################
cat("\n========================================\n")
cat("  (a) Enfoque de agregación ponderada\n")
cat("========================================\n")

# Ajusta los pesos si lo deseas
w_coop <- 0.3
w_trust <- 0.3
w_resou <- 0.2
w_values <- 0.1
w_kinship <- 0.1

agg_matrix <- (m_layers[["Cooperation"]] * w_coop) +
  (m_layers[["Trust"]]       * w_trust) +
  (m_layers[["Resources"]]   * w_resou) +
  (m_layers[["Values"]]      * w_values) +
  (m_layers[["Kinship"]]     * w_kinship)

g_agg <- graph_from_adjacency_matrix(agg_matrix, mode="directed", weighted=TRUE, diag=FALSE)
# Consideramos la versión no dirigida "colapsando" aristas
g_und <- as.undirected(g_agg, mode="collapse", edge.attr.comb="sum")

set.seed(123)  # reproducibilidad
comm_agg <- cluster_louvain(g_und, weights=E(g_und)$weight)
cat("Número de comunidades (agregado):", length(unique(comm_agg$membership)), "\n")
cat("Modularidad (agregado):", modularity(comm_agg), "\n")

df_comm_flat <- data.frame(Nodo=1:vcount(g_und), Comunidad=membership(comm_agg))
cat("\nMuestra de asignaciones de comunidad (10 primeros nodos):\n")
print(head(df_comm_flat, 10))

kbl(head(df_comm_flat,15), format="latex", booktabs=TRUE,
    caption="Asignación de comunidades (agregación ponderada), primeros 15 nodos") %>%
  kable_styling(latex_options="HOLD_position") %>%
  print()

###############################################################################
# (4b) Enfoque especializado con multinet
###############################################################################
cat("\n========================================\n")
cat("  (b) Enfoque especializado con multinet\n")
cat("========================================\n")

omega_val <- 0.2
res_glouvain <- glouvain_ml(net_ml, omega=omega_val)
cat("Número de comunidades con glouvain (omega=", omega_val, "):",
    length(unique(res_glouvain$cid)), "\n")

res_infomap <- infomap_ml(net_ml)
cat("Número de comunidades con infomap (multilayer):",
    length(unique(res_infomap$cid)), "\n")

cat("\n--- Muestra de asignaciones con glouvain_ml() ---\n")
print(head(res_glouvain, 10))

cat("\n--- Muestra de asignaciones con infomap_ml() ---\n")
print(head(res_infomap, 10))

## Tabla de tamaño de comunidades para glouvain
cat("\n*** Resumen de comunidades glouvain_ml() ***\n")
glouvain_summary <- res_glouvain %>%
  group_by(cid) %>%
  summarise(n_vertices = n()) %>%
  arrange(desc(n_vertices))

kbl(glouvain_summary, format="latex", booktabs=TRUE,
    caption="Tamaño de las comunidades con glouvain_ml() en la red multiplex") %>%
  kable_styling(latex_options="HOLD_position") %>%
  print()

###############################################################################
# (4c) Otros métodos de comunidad (ABACUS, Clique Percolation)
###############################################################################
cat("\n========================================\n")
cat("  (c) ABACUS y Clique Percolation\n")
cat("========================================\n")

res_abacus <- abacus_ml(net_ml, min.actors=4, min.layers=2)
cat("\nNúmero de comunidades con ABACUS:", length(unique(res_abacus$cid)), "\n")

res_cclique <- clique_percolation_ml(net_ml, k=4, m=2)
cat("Número de comunidades con Multilayer Clique Percolation:",
    length(unique(res_cclique$cid)), "\n")

###############################################################################
# (5) Tabla comparativa de estadísticas de comunidad (calc_community_stats)
###############################################################################
cat("\n*** Comparando métodos de comunidad con estadísticas descriptivas ***\n")

calc_community_stats <- function(net, comm_df) {
  # Aseguramos que las columnas tengan los nombres correctos
  if(!all(c("actor","layer","cid") %in% names(comm_df))){
    stop("comm_df debe tener columnas actor, layer, cid")
  }
  ncom <- length(unique(comm_df$cid))
  
  all_verts <- nrow(vertices_ml(net))
  assigned_verts <- nrow(dplyr::distinct(comm_df, actor, layer))
  
  all_actors <- num_actors_ml(net)
  assigned_actors <- length(unique(comm_df$actor))
  
  avg_size <- if(ncom>0) assigned_verts / ncom else 0
  
  # (#pares actor–comunidad) / (#actores que están en ≥1 comunidad)
  df_pairs <- dplyr::distinct(comm_df, actor, cid)
  total_pairs <- nrow(df_pairs)
  overlap_level <- if(assigned_actors > 0) total_pairs / assigned_actors else NA
  
  out <- data.frame(
    num = ncom,
    avg_s = avg_size,
    clust_vertices = assigned_verts / all_verts,
    clust_actors = assigned_actors / all_actors,
    actor_overl = overlap_level
  )
  return(out)
}

# Unificar columnas si hiciera falta:
names(res_abacus)  <- c("actor","layer","cid")      # Ajustar según la salida real
names(res_cclique) <- c("actor","layer","cid")      # Ajustar según la salida real
names(res_glouvain) <- c("actor","layer","cid")
names(res_infomap)  <- c("actor","layer","cid")

methods_list <- list(
  abacus  = res_abacus,
  clique  = res_cclique,
  louvain = res_glouvain,
  infomap = res_infomap
)

com_stats <- data.frame()
for(met in names(methods_list)){
  tmp_df <- methods_list[[met]]
  # Verificamos si la tabla está vacía o sin filas
  if(is.null(tmp_df) || nrow(tmp_df) == 0) {
    cat("Aviso: el método", met, "no devolvió comunidades.\n")
    next
  }
  # Calculamos las estadísticas
  tmp_stats <- calc_community_stats(net_ml, tmp_df)
  rownames(tmp_stats) <- met
  com_stats <- rbind(com_stats, tmp_stats)
}

cat("\n--- Tabla de estadísticas de comunidades (ABACUS, CPM, Louvain, Infomap) ---\n")
print(com_stats)

kbl(com_stats, format="latex", booktabs=TRUE,
    caption="Comparación de métodos de comunidad en la red multiplex") %>%
  kable_styling(latex_options="HOLD_position") %>%
  print()






calc_community_stats <- function(net, comm_df, method_name) {
  cat("\nAnalizando estadísticas para método:", method_name, "\n")
  
  # Primero verificamos que tenemos datos válidos
  if(is.null(comm_df) || !is.data.frame(comm_df)) {
    cat("Error: Datos inválidos para", method_name, "\n")
    return(NULL)
  }
  
  # Protegemos cada cálculo en su propio bloque tryCatch
  tryCatch({
    # 1. Estadísticas básicas de la red
    # Obtenemos el número total de vértices y actores de manera robusta
    all_vertices <- vertices_ml(net)
    all_verts <- nrow(all_vertices)
    total_actors <- length(unique(comm_df$actor))
    
    # 2. Número de comunidades
    ncom <- length(unique(comm_df$cid))
    cat("Número de comunidades:", ncom, "\n")
    
    # 3. Vértices asignados (contamos pares actor-layer únicos)
    assigned_vertices <- unique(paste(comm_df$actor, comm_df$layer))
    assigned_verts <- length(assigned_vertices)
    cat("Vértices asignados:", assigned_verts, "\n")
    
    # 4. Actores asignados
    assigned_actors <- length(unique(comm_df$actor))
    cat("Actores asignados:", assigned_actors, "\n")
    
    # 5. Calculamos las métricas finales con protección contra división por cero
    avg_size <- ifelse(ncom > 0, assigned_verts/ncom, 0)
    clust_vertices <- ifelse(all_verts > 0, assigned_verts/all_verts, 0)
    clust_actors <- ifelse(total_actors > 0, assigned_actors/total_actors, 0)
    
    # 6. Calculamos el nivel de superposición
    actor_comm_pairs <- nrow(unique(comm_df[,c("actor", "cid")]))
    actor_overl <- ifelse(assigned_actors > 0, actor_comm_pairs/assigned_actors, 0)
    
    # 7. Creamos el data.frame resultado
    result <- data.frame(
      num = ncom,
      avg_s = avg_size,
      clust_vertices = clust_vertices, 
      clust_actors = clust_actors,
      actor_overl = actor_overl,
      row.names = method_name,
      stringsAsFactors = FALSE
    )
    
    cat("Estadísticas calculadas exitosamente\n")
    return(result)
    
  }, error = function(e) {
    cat("Error en cálculos:", conditionMessage(e), "\n")
    # Retornamos un data.frame con valores por defecto en caso de error
    return(data.frame(
      num = 0,
      avg_s = 0,
      clust_vertices = 0,
      clust_actors = 0,
      actor_overl = 0,
      row.names = method_name,
      stringsAsFactors = FALSE
    ))
  })
}

# Cargar librerías necesarias
suppressMessages({
  library(multinet)
})

# Cargar la red de ejemplo AUCS que viene con el paquete
net <- ml_aucs()

# Ejecutar los diferentes métodos de detección de comunidades
c1 <- abacus_ml(net, 4, 2)
c2 <- clique_percolation_ml(net, 4, 2) 
c3 <- glouvain_ml(net)
c4 <- infomap_ml(net)

# Definir la función corregida para calcular estadísticas
calc_community_stats <- function(net, comm_df, method_name) {
  cat("\nAnalizando estadísticas para método:", method_name, "\n")
  
  # Verificación inicial de datos
  if(is.null(comm_df) || !is.data.frame(comm_df)) {
    cat("Error: Datos inválidos para", method_name, "\n")
    return(NULL)
  }
  
  tryCatch({
    # Estadísticas básicas de la red
    all_vertices <- vertices_ml(net)
    all_verts <- nrow(all_vertices)
    total_actors <- length(unique(comm_df$actor))
    
    # Número de comunidades
    ncom <- length(unique(comm_df$cid))
    cat("Número de comunidades:", ncom, "\n")
    
    # Vértices asignados
    assigned_vertices <- unique(paste(comm_df$actor, comm_df$layer))
    assigned_verts <- length(assigned_vertices)
    cat("Vértices asignados:", assigned_verts, "\n")
    
    # Actores asignados
    assigned_actors <- length(unique(comm_df$actor))
    cat("Actores asignados:", assigned_actors, "\n")
    
    # Métricas finales
    avg_size <- ifelse(ncom > 0, assigned_verts/ncom, 0)
    clust_vertices <- ifelse(all_verts > 0, assigned_verts/all_verts, 0)
    clust_actors <- ifelse(total_actors > 0, assigned_actors/total_actors, 0)
    
    # Nivel de superposición
    actor_comm_pairs <- nrow(unique(comm_df[,c("actor", "cid")]))
    actor_overl <- ifelse(assigned_actors > 0, actor_comm_pairs/assigned_actors, 0)
    
    # Crear data.frame resultado
    result <- data.frame(
      num = ncom,
      avg_s = avg_size,
      clust_vertices = clust_vertices, 
      clust_actors = clust_actors,
      actor_overl = actor_overl,
      row.names = method_name,
      stringsAsFactors = FALSE
    )
    
    cat("Estadísticas calculadas exitosamente\n")
    return(result)
    
  }, error = function(e) {
    cat("Error en cálculos:", conditionMessage(e), "\n")
    return(data.frame(
      num = 0,
      avg_s = 0,
      clust_vertices = 0,
      clust_actors = 0,
      actor_overl = 0,
      row.names = method_name,
      stringsAsFactors = FALSE
    ))
  })
}

# Calcular las estadísticas para cada método
results_list <- list()

# Lista de métodos y sus resultados
methods <- list(
  abacus = c1,
  clique = c2,
  louvain = c3,
  infomap = c4
)

# Procesar cada método
for(met in names(methods)) {
  cat("\n=== Procesando método:", met, "===\n")
  tmp_stats <- calc_community_stats(net, methods[[met]], met)
  if(!is.null(tmp_stats)) {
    results_list[[met]] <- tmp_stats
  }
}

# Combinar resultados
if(length(results_list) > 0) {
  com_stats <- do.call(rbind, results_list)
  cat("\n=== Resumen final de estadísticas ===\n")
  print(com_stats)
} else {
  warning("No se pudieron calcular estadísticas para ningún método.")
}
#######################################################
# (Opcional) Descripción adicional de comunidades
#######################################################
library(dplyr)

average_degree_by_community <- function(net, comm_df) {
  # Calculamos el grado de cada vértice
  vertices_info <- vertices_ml(net)
  colnames(vertices_info) <- c("actor","layer")
  
  deg_values <- degree_ml(net, actors=vertices_info$actor, layers=vertices_info$layer)
  vertices_info$degree <- deg_values
  
  df_joined <- inner_join(comm_df, vertices_info, by=c("actor","layer"))
  
  df_summary <- df_joined %>%
    group_by(cid) %>%
    summarise(n_vertices=n(),
              avg_degree=mean(degree),
              .groups="drop")
  
  return(df_summary)
}

cat("\n--- Ejemplo: estadísticas por comunidad (glouvain) ---\n")
df_g <- average_degree_by_community(net_ml, res_glouvain)
print(df_g)






