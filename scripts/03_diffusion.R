# 0. Librerías necesarias
library(igraph)
library(tidyverse) 
library(ergm)
library(statnet)
library(parallel)
library(patchwork)
library(viridis)

#' Modelo de Difusión de Comportamientos Conflictivos v10
#' Implementa teorías de McAdam (1986), Centola (2007) y Watts (2002)
# improved_diffusion_v10.R

#' Modelo de Difusión de Comportamientos Conflictivos
#' Implementa teorías de McAdam (1986), Centola (2007) y Watts (2002)
# Versión mejorada de la función principal
improved_diffusion_v10 <- function(networks, attributes, seeds,
                                   params = default_params,
                                   max_iterations = 1000,
                                   convergence_threshold = 0.001) {
  
  # Inicialización
  n_nodes <- nrow(networks[[1]])
  states <- rep(0, n_nodes)
  states[seeds] <- 1
  
  # Matrices de memoria
  history <- matrix(0, nrow = max_iterations, ncol = n_nodes)
  history[1,] <- states
  
  exposure_memory <- matrix(0, nrow = max_iterations, ncol = n_nodes)
  clustering_memory <- matrix(0, nrow = max_iterations, ncol = n_nodes)
  
  # Proceso de difusión mejorado
  for(iter in 2:max_iterations) {
    old_states <- states
    
    # Procesar nodos inactivos
    inactive_nodes <- which(states == 0)
    
    for(i in inactive_nodes) {
      # Calcular influencia con el nuevo método
      influence <- calc_complex_influence(i, networks, states, params$weights)
      
      # Calcular efectos de memoria y clustering
      memory_effect <- mean(exposure_memory[1:max(1,iter-1), i], na.rm=TRUE)
      cluster_effect <- calc_local_clustering(i, networks, params$weights)
      
      # Calcular umbral dinámico
      threshold <- get_dynamic_threshold(i, attributes, params, memory_effect, cluster_effect)
      
      # Regla de adopción
      if(influence >= threshold) {
        adoption_prob <- params$adoption_friction * (1 + memory_effect)
        states[i] <- rbinom(1, 1, min(1, adoption_prob))
      }
    }
    
    # Actualizar memorias
    history[iter,] <- states
    exposure_memory[iter,] <- update_exposure_memory(states, networks)
    clustering_memory[iter,] <- sapply(1:n_nodes, 
                                       function(i) calc_local_clustering(i, networks, params$weights))
    
    # Verificar convergencia
    if(check_convergence(history, iter, convergence_threshold)) break
  }
  
  return(process_results(states, history, iter, attributes, networks, params))
}


# auxiliary_functions.R

#' Cálculo de influencia compleja siguiendo a Centola & Macy (2007)
#' Implementa el concepto de contagio complejo donde múltiples exposiciones 
#' son necesarias para la adopción
# Función mejorada para calcular influencia compleja
calc_complex_influence <- function(node, networks, states, weights) {
  influence <- 0
  total_possible <- 0
  
  # Calcular influencia ponderada por tipo de vínculo
  for(n in seq_along(networks)) {
    neighbors <- which(networks[[n]][node,] > 0)
    if(length(neighbors) > 0) {
      # Ponderación por tipo de vínculo (siguiendo a Diani)
      active_neighbors <- sum(states[neighbors])
      influence <- influence + (active_neighbors * weights[n])
      total_possible <- total_possible + (length(neighbors) * weights[n])
    }
  }
  
  # Normalizar y aplicar factor de contagio complejo
  influence_norm <- if(total_possible > 0) influence/total_possible else 0
  
  # Aplicar umbral de redundancia según Centola
  return(ifelse(influence_norm >= 0.2, influence_norm, 0))
}

# Función mejorada para calcular clustering local
calc_local_clustering <- function(node, networks, weights) {
  cluster_scores <- sapply(networks, function(x) {
    neighbors <- which(x[node,] > 0)
    if(length(neighbors) < 2) return(0)
    submat <- x[neighbors, neighbors]
    clustering <- sum(submat) / (length(neighbors) * (length(neighbors)-1))
    return(clustering * weights[match(x, networks)])
  })
  return(mean(cluster_scores))
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
#' Implementa teoría de Granovetter sobre umbrales heterogéneos
#' Función get_threshold corregida
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
                                     density_range = seq(0.01, 0.2, by = 0.01)) {
  results <- list()
  
  for(d in density_range) {
    # Crear redes con densidad d manteniendo estructura
    sparse_networks <- lapply(base_networks, function(net) {
      threshold <- quantile(net[net > 0], 1 - d)
      net * (net >= threshold)
    })
    
    # Correr simulación
    sim_result <- improved_diffusion_v10(sparse_networks, attributes, seeds)
    results[[as.character(d)]] <- sim_result
  }
  
  return(results)
}

#' Función para variar fracción de umbrales bajos
test_threshold_sensitivity <- function(base_thresholds,
                                       low_threshold_props = seq(0.1, 0.5, by = 0.05)) {
  results <- list()
  
  for(p in low_threshold_props) {
    # Modificar distribución de umbrales
    modified_thresholds <- base_thresholds
    n_low <- round(length(base_thresholds) * p)
    modified_thresholds[1:n_low] <- modified_thresholds[1:n_low] * 0.5
    
    # Correr simulación
    sim_result <- improved_diffusion_v10(networks, attributes, seeds,
                                         thresholds = modified_thresholds)
    results[[as.character(p)]] <- sim_result
  }
  
  return(results)
}

#' Función para variar nivel de represión
test_repression_sensitivity <- function(repression_levels = seq(0, 0.5, by = 0.05)) {
  results <- list()
  
  for(r in repression_levels) {
    sim_result <- improved_diffusion_v10(networks, attributes, seeds,
                                         repression_level = r)
    results[[as.character(r)]] <- sim_result
  }
  
  return(results)
}


# Parámetros del modelo ajustados
params <- list(
  weights = c(0.35, 0.35, 0.15, 0.15), # Mayor peso a confianza y valores
  thresholds = list(
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
  # Nuevos parámetros de difusión
  adoption_friction = 0.3,     # Fricción en adopción
  memory_decay = 0.05,        # Decay de memoria
  clustering_impact = 0.2,     # Impacto de clustering
  repression_base = 0.1,      # Represión base
  identity_bonus = 0.15       # Bonus por identidad compartida
)

# Script Principal de Análisis
# Basado en teorías de movimientos sociales y contagio complejo





# 1. Preparación del ambiente ----
library(tidyverse)
library(igraph)
library(ergm)
library(statnet)
library(parallel)
library(patchwork)
library(viridis)

# Cargar funciones auxiliares
source("blog/posts/07-social_movements/improved_diffusion_v10.R")
source("blog/posts/07-social_movements/auxiliary_functions.R")

# 2. Carga y Preparación de Datos ----
# Cargar redes
coop_matrix      <- as.matrix(read.csv("blog/posts/07-social_movements/datos/CoopNet.csv", header=FALSE))
trust_matrix     <- as.matrix(read.csv("blog/posts/07-social_movements/datos/ConfianzaNet.csv", header=FALSE))
resources_matrix <- as.matrix(read.csv("blog/posts/07-social_movements/datos/RecursosNet.csv", header=FALSE))
values_matrix    <- as.matrix(read.csv("blog/posts/07-social_movements/datos/ValoresNet.csv", header=FALSE))

# Cargar atributos
org_attributes <- read.csv("blog/posts/07-social_movements/datos/Atributos _org2011.csv", header=TRUE) %>%
  mutate(
    conflictivo = trimws(Orientación) == "Conflictiva",
    tipo = Tipo_de_organización,
    ubicacion = Ubicación,
    orientacion = Orientación,
    block = as.factor(X1PosciónCONCOR)
  )

# Crear lista de redes
networks <- list(
  trust = trust_matrix,    # Red de confianza (peso 0.35)
  values = values_matrix,   # Red de valores (peso 0.35)
  coop = coop_matrix,     # Red de cooperación (peso 0.15)
  resources = resources_matrix # Red de recursos (peso 0.15)
)

# 5. Simulaciones de Difusión ----
# Definir semillas (organizaciones políticas de pobladores)
seeds <- which(org_attributes$tipo == " Organización política de pobladores ")

# Simulación base
base_results <- improved_diffusion_v10(
  networks = networks,
  attributes = org_attributes,
  seeds = seeds,
  weights = c(0.35, 0.35, 0.15, 0.15),
  deactivation_threshold = 0.3,
  repression_level = 0.1
)

# 6. Análisis de Sensibilidad ----
# Variar densidad de red
density_results <- test_density_sensitivity(
  base_networks = networks, 
  attributes = org_attributes, 
  seeds = seeds
)

# Variar umbrales
threshold_results <- test_threshold_sensitivity(
  networks = networks,
  attributes = org_attributes, 
  seeds = seeds,
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
)


# Variar represión
repression_results <- test_repression_sensitivity(
  networks = networks,
  attributes = org_attributes, 
  seeds = seeds
)

# 7. Visualización de Resultados ----
# Función para visualizar trayectorias de difusión
plot_diffusion_trajectories <- function(results) {
  df_trajectories <- data.frame(
    iteration = 1:results$n_iterations,
    adoption = results$final_metrics$change_trajectory
  )
  
  ggplot(df_trajectories, aes(x = iteration, y = adoption)) +
    geom_line(size = 1.2, color = "#2C3E50") +
    geom_point(size = 3, color = "#E74C3C") +
    theme_minimal() +
    labs(
      title = "Trayectoria de Difusión",
      x = "Iteración",
      y = "Proporción de Adopción"
    )
}

# Visualizar resultados de sensibilidad
plot_sensitivity_analysis <- function(density_results, 
                                      threshold_results,
                                      repression_results) {
  # Crear dataframes para cada análisis
  df_density <- map_df(density_results, ~.x$final_metrics$adoption_rate,
                       .id = "density")
  df_threshold <- map_df(threshold_results, ~.x$final_metrics$adoption_rate,
                         .id = "threshold_prop")
  df_repression <- map_df(repression_results, ~.x$final_metrics$adoption_rate,
                          .id = "repression")
  
  # Crear visualizaciones
  p1 <- ggplot(df_density, aes(x = as.numeric(density), 
                               y = adoption_rate)) +
    geom_line() +
    labs(title = "Efecto de Densidad de Red")
  
  p2 <- ggplot(df_threshold, aes(x = as.numeric(threshold_prop),
                                 y = adoption_rate)) +
    geom_line() +
    labs(title = "Efecto de Proporción de Umbrales Bajos")
  
  p3 <- ggplot(df_repression, aes(x = as.numeric(repression),
                                  y = adoption_rate)) +
    geom_line() +
    labs(title = "Efecto de Nivel de Represión")
  
  p1 + p2 + p3 + plot_layout(ncol = 1)
}

# 8. Guardar Resultados ----
saveRDS(list(
  network_metrics = network_metrics,
  ergm_model = ergm_model,
  base_results = base_results,
  sensitivity = list(
    density = density_results,
    threshold = threshold_results,
    repression = repression_results
  )
), "resultados_completos.rds")

# 9. Generar Reporte ----
# Ver archivo separate_report_generator.R






# 7. Visualización de Resultados ----

# Función para visualizar trayectorias de difusión
plot_diffusion_trajectories <- function(results) {
  # Verificar que el objeto de resultados tenga la estructura esperada
  if (!"final_metrics" %in% names(results) || 
      !"change_trajectory" %in% names(results$final_metrics)) {
    stop("Los resultados no tienen la estructura esperada. Verifique la simulación.")
  }
  
  # Preparar datos de trayectoria
  df_trajectories <- data.frame(
    iteration = 1:results$n_iterations,
    adoption = results$final_metrics$change_trajectory
  )
  
  # Crear el gráfico con ggplot2
  ggplot(df_trajectories, aes(x = iteration, y = adoption)) +
    geom_line(size = 1.2, color = "#2C3E50") +
    geom_point(size = 3, color = "#E74C3C") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    labs(
      title = "Trayectoria de Difusión de Comportamientos",
      x = "Número de Iteraciones",
      y = "Proporción de Adopción"
    )
}

# Función para visualizar análisis de sensibilidad
plot_sensitivity_analysis <- function(density_results, 
                                      threshold_results,
                                      repression_results) {
  # Función auxiliar para extraer tasas de adopción
  extract_adoption_rates <- function(results, id_name) {
    # Manejar diferentes estructuras de resultados
    if ("summary" %in% names(results)) {
      # Si tiene estructura de test_threshold_sensitivity
      df <- results$summary
      names(df)[names(df) == id_name] <- "x"
      names(df)[names(df) == "adoption_rate"] <- "y"
    } else {
      # Estructura de otros análisis de sensibilidad
      df <- map_df(results, ~.x$final_metrics$adoption_rate, .id = id_name)
      names(df)[names(df) == id_name] <- "x"
      names(df)[names(df) == "final_metrics.adoption_rate"] <- "y"
    }
    
    # Convertir x a numérico
    df$x <- as.numeric(as.character(df$x))
    
    return(df)
  }
  
  # Extraer datos
  df_density <- extract_adoption_rates(density_results, "density")
  df_threshold <- extract_adoption_rates(threshold_results, "threshold_prop")
  df_repression <- extract_adoption_rates(repression_results, "repression")
  
  # Crear visualizaciones
  p1 <- ggplot(df_density, aes(x = x, y = y)) +
    geom_line(color = "#3498DB", size = 1.2) +
    geom_point(color = "#2980B9", size = 3) +
    theme_minimal() +
    labs(
      title = "Efecto de Densidad de Red",
      x = "Densidad de Red",
      y = "Tasa de Adopción"
    )
  
  p2 <- ggplot(df_threshold, aes(x = x, y = y)) +
    geom_line(color = "#2ECC71", size = 1.2) +
    geom_point(color = "#27AE60", size = 3) +
    theme_minimal() +
    labs(
      title = "Efecto de Proporción de Umbrales Bajos",
      x = "Proporción de Reducción de Umbral",
      y = "Tasa de Adopción"
    )
  
  p3 <- ggplot(df_repression, aes(x = x, y = y)) +
    geom_line(color = "#E74C3C", size = 1.2) +
    geom_point(color = "#C0392B", size = 3) +
    theme_minimal() +
    labs(
      title = "Efecto de Nivel de Represión",
      x = "Nivel de Represión",
      y = "Tasa de Adopción"
    )
  
  # Combinar gráficos
  library(patchwork)
  return(p1 / p2 / p3)  # Usar layout vertical
}

# Ejemplo de uso
# Asumiendo que ya has corrido las simulaciones
# plot_diffusion_trajectories(base_results)
# plot_sensitivity_analysis(
#   density_results, 
#   threshold_results, 
#   repression_results
# )


# Ejemplo de uso
# Asumiendo que ya has corrido las simulaciones
plot_diffusion_trajectories(base_results)
 plot_sensitivity_analysis(
   density_results, 
   threshold_results, 
   repression_results
 )

 
 
 # Check the structure of your results
 str(density_results)
 str(threshold_results)
 str(repression_results)

 
 
 # Diagnostic function to understand the structure of results
 diagnose_sensitivity_results <- function(density_results, 
                                          threshold_results, 
                                          repression_results) {
   cat("Density Results Structure:\n")
   print(str(density_results))
   
   cat("\nThreshold Results Structure:\n")
   print(str(threshold_results))
   
   cat("\nRepression Results Structure:\n")
   print(str(repression_results))
 }
 
 # Modified extraction function with error handling
 extract_adoption_rates <- function(results) {
   # Check if results is empty or not a list
   if (length(results) == 0 || !is.list(results)) {
     return(data.frame(x = numeric(), y = numeric()))
   }
   
   # If results have a 'summary' component (from test_threshold_sensitivity)
   if ("summary" %in% names(results)) {
     return(results$summary)
   }
   
   # For other result types
   adoption_rates <- lapply(results, function(result) {
     # Safely extract adoption rate, return NA if not found
     tryCatch(
       result$final_metrics$adoption_rate, 
       error = function(e) NA
     )
   })
   
   # Remove NA values
   adoption_rates <- adoption_rates[!is.na(adoption_rates)]
   
   # Create data frame
   data.frame(
     x = as.numeric(names(adoption_rates)),
     y = unlist(adoption_rates)
   )
 }
 
 # Robust visualization function
 plot_sensitivity_analysis <- function(density_results, 
                                       threshold_results,
                                       repression_results) {
   # Extract rates with error handling
   df_density <- extract_adoption_rates(density_results)
   df_threshold <- extract_adoption_rates(threshold_results)
   df_repression <- extract_adoption_rates(repression_results)
   
   # Debugging output
   cat("Density DataFrame:\n")
   print(df_density)
   cat("\nThreshold DataFrame:\n")
   print(df_threshold)
   cat("\nRepression DataFrame:\n")
   print(df_repression)
   
   # Require ggplot2 and patchwork
   library(ggplot2)
   library(patchwork)
   
   # Plots with error handling
   plots <- list()
   
   # Density plot
   if (nrow(df_density) > 0) {
     plots$density <- ggplot(df_density, aes(x = x, y = y)) +
       geom_line(color = "#3498DB", size = 1.2) +
       geom_point(color = "#2980B9", size = 3) +
       theme_minimal() +
       labs(
         title = "Efecto de Densidad de Red",
         x = "Densidad de Red",
         y = "Tasa de Adopción"
       )
   }
   
   # Threshold plot
   if (nrow(df_threshold) > 0) {
     plots$threshold <- ggplot(df_threshold, aes(x = x, y = y)) +
       geom_line(color = "#2ECC71", size = 1.2) +
       geom_point(color = "#27AE60", size = 3) +
       theme_minimal() +
       labs(
         title = "Efecto de Proporción de Umbrales Bajos",
         x = "Proporción de Reducción de Umbral",
         y = "Tasa de Adopción"
       )
   }
   
   # Repression plot
   if (nrow(df_repression) > 0) {
     plots$repression <- ggplot(df_repression, aes(x = x, y = y)) +
       geom_line(color = "#E74C3C", size = 1.2) +
       geom_point(color = "#C0392B", size = 3) +
       theme_minimal() +
       labs(
         title = "Efecto de Nivel de Represión",
         x = "Nivel de Represión",
         y = "Tasa de Adopción"
       )
   }
   
   # Combine available plots
   if (length(plots) > 0) {
     # Use patchwork to combine plots vertically
     do.call(c, plots)
   } else {
     stop("No valid plots could be created. Please check your sensitivity analysis results.")
   }
 }
 
 # First, diagnose the results
 diagnose_sensitivity_results(
   density_results, 
   threshold_results, 
   repression_results
 )
 
 # Then attempt to plot
 plot_sensitivity_analysis(
   density_results, 
   threshold_results, 
   repression_results
 )
 
 
 
 

 
 # Visualización de análisis de sensibilidad
 plot_sensitivity_analysis <- function(density_results, 
                                       threshold_results,
                                       repression_results) {
   library(ggplot2)
   library(patchwork)
   library(dplyr)
   
   # Preparar datos de densidad
   density_data <- data.frame(
     density = as.numeric(names(density_results)),
     adoption_rate = sapply(density_results, function(res) res$final_metrics$adoption_rate)
   )
   
   # Gráfico de Densidad de Red
   density_plot <- ggplot(density_data, aes(x = density, y = adoption_rate)) +
     geom_line(color = "#3498DB", size = 1.2) +
     geom_point(color = "#2980B9", size = 3) +
     theme_minimal() +
     labs(
       title = "Efecto de Densidad de Red",
       x = "Densidad de Red",
       y = "Tasa de Adopción"
     )
   
   # Gráfico de Umbrales (usando summary de threshold_results)
   threshold_plot <- ggplot(threshold_results$summary, 
                            aes(x = threshold_reduction, y = adoption_rate)) +
     geom_line(color = "#2ECC71", size = 1.2) +
     geom_point(color = "#27AE60", size = 3) +
     theme_minimal() +
     labs(
       title = "Efecto de Proporción de Umbrales Bajos",
       x = "Proporción de Reducción de Umbral",
       y = "Tasa de Adopción"
     )
   
   # Preparar datos de represión
   repression_data <- data.frame(
     repression = as.numeric(names(repression_results)),
     adoption_rate = sapply(repression_results, function(res) res$final_metrics$adoption_rate)
   )
   
   # Gráfico de Represión
   repression_plot <- ggplot(repression_data, 
                             aes(x = repression, y = adoption_rate)) +
     geom_line(color = "#E74C3C", size = 1.2) +
     geom_point(color = "#C0392B", size = 3) +
     theme_minimal() +
     labs(
       title = "Efecto de Nivel de Represión",
       x = "Nivel de Represión",
       y = "Tasa de Adopción"
     )
   
   # Combinar plots
   combined_plot <- density_plot / threshold_plot / repression_plot
   
   # Imprimir plots
   print(combined_plot)
   
   return(combined_plot)
 }
 
 # Ejecutar la función de visualización
 sensitivity_plots <- plot_sensitivity_analysis(
   density_results, 
   threshold_results, 
   repression_results
 )
 