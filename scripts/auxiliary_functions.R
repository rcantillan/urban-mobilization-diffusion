
# auxiliary_functions.R

#' Cálculo de influencia compleja siguiendo a Centola & Macy (2007)
#' Implementa el concepto de contagio complejo donde múltiples exposiciones 
#' son necesarias para la adopción
calc_complex_influence <- function(node, networks, states, weights) {
  influence <- 0
  total_possible <- 0
  
  # Calcular influencia ponderada por tipo de vínculo
  for(n in seq_along(networks)) {
    neighbors <- which(networks[[n]][node,] > 0)
    if(length(neighbors) > 0) {
      # Ponderación por tipo de vínculo siguiendo a Diani
      tie_weight <- weights[n]
      
      # Suma ponderada de estados de vecinos
      active_neighbors <- sum(states[neighbors])
      influence <- influence + (active_neighbors * tie_weight)
      total_possible <- total_possible + (length(neighbors) * tie_weight)
    }
  }
  
  # Normalizar influencia total
  return(if(total_possible > 0) influence/total_possible else 0)
}

#' Cálculo de presión de desactivación
#' Implementa mecanismo de salida por represión o aislamiento
#' Basado en McAdam (1986) sobre costos de activismo de alto riesgo
calc_deactivation_pressure <- function(node, networks, states, repression_level) {
  # Presión base por represión
  pressure <- repression_level
  
  # Añadir presión por aislamiento
  active_neighbors <- 0
  total_neighbors <- 0
  
  for(net in networks) {
    neighbors <- which(net[node,] > 0)
    if(length(neighbors) > 0) {
      active_neighbors <- active_neighbors + sum(states[neighbors])
      total_neighbors <- total_neighbors + length(neighbors)
    }
  }
  
  # Calcular ratio de aislamiento
  isolation_ratio <- if(total_neighbors > 0) {
    1 - (active_neighbors/total_neighbors)
  } else {
    1
  }
  
  return(pressure * (1 + isolation_ratio))
}

#' Cálculo de umbral dinámico
get_threshold <- function(node, attributes, base_thresholds, exposure_effect, cluster_effect,
                          use_heterogeneous_thresholds = TRUE) {
  
  if(use_heterogeneous_thresholds) {
    base_threshold <- base_thresholds$by_type[attributes$tipo[node]]
  } else {
    base_threshold <- base_thresholds$base
  }
  
  # Ajuste dinámico del umbral
  return(base_threshold * 
           exp(-0.05 * exposure_effect) * 
           (1 - 0.2 * cluster_effect))
}


#' Actualización de memoria de clustering
#' Calcula el efecto de clustering en redes multiplex
#' Implementa la idea de que la adopción depende de la densidad del vecindario
update_clustering_memory <- function(states, networks, weights) {
  sapply(1:length(states), function(node) {
    # Calcular clustering por cada capa de red
    layer_clusterings <- sapply(seq_along(networks), function(n) {
      # Obtener vecinos del nodo en esta capa
      neighbors <- which(networks[[n]][node,] > 0)
      
      if(length(neighbors) > 0) {
        # Calcular proporción de vecinos activos
        active_neighbors <- sum(states[neighbors])
        clustering_coef <- active_neighbors / length(neighbors)
        
        # Ponderar por el peso de la capa
        return(clustering_coef * weights[n])
      } else {
        return(0)
      }
    })
    
    # Combinar clustering de todas las capas
    return(mean(layer_clusterings, na.rm = TRUE))
  })
}


#' Actualización de memoria de exposición
#' Implementa teoría de memoria colectiva en movimientos sociales
update_exposure_memory <- function(states, node_neighbors) {
  sapply(seq_along(states), function(i) {
    if(length(node_neighbors[[i]]) > 0) {
      mean(states[node_neighbors[[i]]])
    } else {
      0
    }
  })
}

#' Verificación de convergencia
#' Basado en teoría de Watts sobre cascadas globales
check_convergence <- function(history, iter, threshold) {
  if(iter > 20) {
    recent_window <- (iter-19):iter
    change_rate <- mean(diff(colMeans(history[recent_window,])))
    return(abs(change_rate) < threshold)
  }
  return(FALSE)
}

#' Cálculo de métricas finales
#' Implementa indicadores clave según literatura de movimientos sociales

calculate_final_metrics <- function(states, history, iter, attributes = NULL, role_scores = NULL) {
  list(
    # Tamaño final de adopción (Watts)
    adoption_rate = mean(states),
    
    # Velocidad de difusión (Centola)
    time_to_equilibrium = iter,
    
    # Trayectoria de cambio (McAdam)
    change_trajectory = rowMeans(history[1:iter,]),
    
    # Patrones por tipo y rol (solo si attributes y role_scores están presentes)
    adoption_by_type = if(!is.null(attributes)) {
      tapply(states, attributes$tipo, mean)
    } else {
      NULL
    },
    adoption_by_role = if(!is.null(role_scores)) {
      tapply(states, role_scores$role, mean)
    } else {
      NULL
    }
  )
}



#' Análisis de Sensibilidad
#' Explora sistemáticamente el espacio de parámetros siguiendo a Watts (2002)

#' Función para variar densidad de red
test_density_sensitivity <- function(base_networks, 
                                     attributes,  # Add this parameter
                                     seeds,       # Add this parameter
                                     density_range = seq(0.01, 0.2, by = 0.01)) {
  results <- list()
  
  for(d in density_range) {
    # Crear redes con densidad d manteniendo estructura
    sparse_networks <- lapply(base_networks, function(net) {
      threshold <- quantile(net[net > 0], 1 - d)
      net * (net >= threshold)
    })
    
    # Correr simulación
    sim_result <- improved_diffusion_v10(
      networks = sparse_networks, 
      attributes = attributes, 
      seeds = seeds
    )
    results[[as.character(d)]] <- sim_result
  }
  
  return(results)
}


#' Análisis de Sensibilidad de Umbrales
#' 
#' Esta función explora cómo diferentes configuraciones de umbrales 
#' afectan la difusión de comportamientos en redes sociales
#'
#' @param networks Lista de matrices de redes multiplex
#' @param attributes Dataframe con atributos de los nodos
#' @param seeds Índices de nodos semilla que inician la difusión
#' @param base_thresholds Umbrales base para la adopción de comportamientos
#' @param low_threshold_props Proporción de reducción de umbrales a probar
#'
#' @return Lista de resultados de simulaciones con diferentes configuraciones de umbrales

test_threshold_sensitivity <- function(
    networks,           # Redes multiplex de entrada
    attributes,         # Atributos de los nodos
    seeds,              # Nodos semilla
    base_thresholds = list(  # Umbrales por defecto
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
    ),
    low_threshold_props = seq(0.1, 0.5, by = 0.05)  # Proporciones de reducción de umbral
) {
  # Contenedor para resultados de simulaciones
  results <- list()
  
  # Iterar sobre diferentes proporciones de reducción de umbral
  for(p in low_threshold_props) {
    # Crear una copia modificable de los umbrales base
    modified_thresholds <- base_thresholds
    
    # Modificar umbrales por tipo de organización
    # Reducir los umbrales proporcionalmente a 'p'
    modified_thresholds$by_type <- sapply(
      modified_thresholds$by_type, 
      function(x) x * (1 - p)  # Reducir umbrales
    )
    
    # Modificar también el umbral base
    modified_thresholds$base <- base_thresholds$base * (1 - p)
    
    # Ejecutar simulación con umbrales modificados
    sim_result <- improved_diffusion_v10(
      networks = networks, 
      attributes = attributes, 
      seeds = seeds,
      base_thresholds = modified_thresholds
    )
    
    # Almacenar resultados con la proporción de reducción como identificador
    results[[as.character(p)]] <- list(
      # Métricas principales
      adoption_rate = sim_result$final_metrics$adoption_rate,
      time_to_equilibrium = sim_result$final_metrics$time_to_equilibrium,
      
      # Adopción por tipo de organización
      adoption_by_type = sim_result$final_metrics$adoption_by_type,
      
      # Parámetros utilizados
      threshold_reduction = p,
      modified_thresholds = modified_thresholds,
      
      # Resultado completo de la simulación
      full_simulation = sim_result
    )
  }
  
  # Preparar datos para visualización
  results_summary <- data.frame(
    threshold_reduction = as.numeric(names(results)),
    adoption_rate = sapply(results, `[[`, "adoption_rate"),
    time_to_equilibrium = sapply(results, `[[`, "time_to_equilibrium")
  )
  
  # Retornar resultados completos y resumen
  return(list(
    detailed_results = results,
    summary = results_summary
  ))
}

#' Función para variar nivel de represión
test_repression_sensitivity <- function(networks,    # Add these parameters
                                        attributes,  # 
                                        seeds,       # 
                                        repression_levels = seq(0, 0.5, by = 0.05)) {
  results <- list()
  
  for(r in repression_levels) {
    sim_result <- improved_diffusion_v10(
      networks = networks, 
      attributes = attributes, 
      seeds = seeds,
      repression_level = r
    )
    results[[as.character(r)]] <- sim_result
  }
  
  return(results)
}
