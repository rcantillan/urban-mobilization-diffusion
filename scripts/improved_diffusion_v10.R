#' Modelo de Difusión de Comportamientos Conflictivos v10
#' Implementa teorías de McAdam (1986), Centola (2007) y Watts (2002)
#' Modelo de Difusión de Comportamientos Conflictivos v10
#' 
#' Esta función simula la difusión de comportamientos en redes multiplex, 
#' integrando teorías de McAdam (1986), Centola (2007) y Watts (2002)
#'
#' @param networks Lista de matrices de redes multiplex (confianza, valores, cooperación, recursos)
#' @param attributes Dataframe con atributos de los nodos
#' @param seeds Índices de nodos semilla que inician la difusión
#' @param weights Pesos para cada capa de red (confianza, valores, cooperación, recursos)
#' @param deactivation_threshold Umbral de desactivación para nodos activos
#' @param repression_level Nivel de represión que afecta la permanencia de nodos activos
#' @param max_iterations Número máximo de iteraciones de simulación
#' @param convergence_threshold Umbral para determinar convergencia del modelo
#' @param memory_window Ventana de memoria para efectos de exposición y clustering
#' @param use_heterogeneous_thresholds Usar umbrales diferenciados por tipo de organización
#' @param base_thresholds Umbrales base para la adopción de comportamientos
#'
#' @return Lista con resultados de la simulación: estados finales, historia, métricas, etc.
improved_diffusion_v10 <- function(
    networks, 
    attributes, 
    seeds,
    weights = c(0.35, 0.35, 0.15, 0.15),
    deactivation_threshold = 0.3,
    repression_level = 0.1,
    max_iterations = 1000,
    convergence_threshold = 0.001,
    memory_window = 5,
    use_heterogeneous_thresholds = TRUE,
    base_thresholds = list(
      base = 0.05,
      by_type = c(
        " Club deportivo " = 0.15,
        " Comité de vivienda " = 0.06,
        " Organización cultural " = 0.10,
        " Organización política de base " = 0.07,
        " Organización política de pobladores " = 0.05,
        " Organización vecinal " = 0.10,
        "Otros" = 0.15
      )
    )
) {
  
  # 1. Validación y preprocesamiento de redes ----
  # Convertir todas las redes a matrices
  networks <- lapply(networks, as.matrix)
  n_nodes <- nrow(networks[[1]])
  
  # 2. Crear red compuesta ponderada ----
  # Siguiendo el enfoque de Diani sobre redes multiplex
  composite_net <- Reduce("+", Map("*", networks, weights))
  
  # 3. Calcular medidas estructurales ----
  g <- igraph::graph_from_adjacency_matrix(composite_net, 
                                           weighted = TRUE,
                                           mode = "directed")
  
  # Calcular medidas de centralidad para identificación de roles
  centrality_measures <- list(
    between = igraph::betweenness(g, normalized = TRUE),
    degree = igraph::degree(g, mode = "all", normalized = TRUE),
    eigen = igraph::eigen_centrality(g)$vector,
    authority = igraph::authority_score(g)$vector
  )
  
  # 4. Identificar roles estructurales ----
  role_scores <- data.frame(
    node = 1:n_nodes,
    between_score = scale(centrality_measures$between),
    degree_score = scale(centrality_measures$degree),
    eigen_score = scale(centrality_measures$eigen),
    authority_score = scale(centrality_measures$authority)
  )
  
  # Calcular scores compuestos para roles
  role_scores$broker_score <- role_scores$between_score + role_scores$degree_score
  role_scores$core_score <- role_scores$eigen_score + role_scores$authority_score
  
  # Asignar roles según percentiles
  role_scores$role <- with(role_scores, {
    ifelse(broker_score > quantile(broker_score, 0.85), "broker",
           ifelse(core_score > quantile(core_score, 0.85), "core",
                  ifelse(authority_score > quantile(authority_score, 0.85), 
                         "authority", "peripheral")))
  })
  
  # 5. Crear lista de vecinos para cálculos de influencia ----
  node_neighbors <- lapply(1:n_nodes, function(i) {
    unique(unlist(lapply(networks, function(x) which(x[i,] > 0))))
  })
  
  # 6. Inicialización de estados y memorias ----
  states <- rep(0, n_nodes)
  states[seeds] <- 1
  
  # Matrices para rastrear historia y memoria
  history <- matrix(0, nrow = max_iterations, ncol = n_nodes)
  history[1,] <- states
  
  exposure_memory <- matrix(0, nrow = max_iterations, ncol = n_nodes)
  exposure_memory[1,] <- update_exposure_memory(states, node_neighbors)
  
  clustering_memory <- matrix(0, nrow = max_iterations, ncol = n_nodes)
  clustering_memory[1,] <- update_clustering_memory(states, networks, weights)
  
  # 7. Proceso de difusión ----
  for(iter in 2:max_iterations) {
    old_states <- states
    
    # Procesar todos los nodos en cada iteración
    for(i in 1:n_nodes) {
      # Calcular influencia compleja según Centola
      influence <- calc_complex_influence(i, networks, states, weights)
      
      # Calcular efectos para umbral dinámico
      exposure_effect <- mean(exposure_memory[1:max(1,iter-1), i], na.rm=TRUE)
      cluster_effect <- mean(clustering_memory[1:max(1,iter-1), i], na.rm=TRUE)
      
      # Obtener umbral dinámico
      threshold <- get_threshold(
        node = i,
        attributes = attributes,
        base_thresholds = base_thresholds,
        exposure_effect = exposure_effect,
        cluster_effect = cluster_effect,
        use_heterogeneous_thresholds = use_heterogeneous_thresholds
      )
      
      # 8. Reglas de activación/desactivación ----
      # Nodo inactivo puede activarse
      if(states[i] == 0 && influence >= threshold) {
        # Probabilidad de adopción con fricción
        adoption_prob <- 0.3 * (1 + exposure_effect)
        states[i] <- rbinom(1, 1, min(1, adoption_prob))
      } 
      # Nodo activo puede desactivarse
      else if(states[i] == 1) {
        # Verificar posible desactivación
        deactivation_pressure <- calc_deactivation_pressure(
          i, networks, states, repression_level
        )
        if(deactivation_pressure > deactivation_threshold) {
          states[i] <- 0
        }
      }
    }
    
    # 9. Actualizar memorias del proceso ----
    history[iter,] <- states
    exposure_memory[iter,] <- update_exposure_memory(states, node_neighbors)
    clustering_memory[iter,] <- update_clustering_memory(states, networks, weights)
    
    # 10. Verificar convergencia ----
    if(check_convergence(history, iter, convergence_threshold)) break
  }
  
  # 11. Calcular métricas finales ----
  final_metrics <- calculate_final_metrics(
    states, history, iter, attributes, role_scores
  )
  
  # 12. Retornar resultados completos ----
  return(list(
    final_states = states,
    history = history[1:iter,],
    n_iterations = iter,
    role_scores = role_scores,
    centrality_measures = centrality_measures,
    exposure_memory = exposure_memory[1:iter,],
    clustering_memory = clustering_memory[1:iter,],
    final_metrics = final_metrics,
    parameters = list(
      thresholds = base_thresholds,
      weights = weights,
      repression_level = repression_level,
      deactivation_threshold = deactivation_threshold,
      use_heterogeneous_thresholds = use_heterogeneous_thresholds
    )
  ))
}

