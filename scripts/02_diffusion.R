
# Cargar librerías necesarias
library(tidyverse)
library(igraph)

# Cargar las matrices de redes desde los archivos CSV
conf_net <- as.matrix(read.csv("blog/posts/07-social_movements/datos/ConfianzaNet.csv", header=FALSE))
coop_net <- as.matrix(read.csv("blog/posts/07-social_movements/datos/CoopNet.csv", header=FALSE))
rec_net <-  as.matrix(read.csv("blog/posts/07-social_movements/datos/RecursosNet.csv", header=FALSE))
val_net <-  as.matrix(read.csv("blog/posts/07-social_movements/datos/ValoresNet.csv", header=FALSE))

# Cargar atributos
attr <- read.csv("blog/posts/07-social_movements/datos/Atributos _org2011.csv") %>%
  mutate(
    bloque_concor = X1PosciónCONCOR,
    tipo = Tipo,
    conflictivo = trimws(Orientación) == "Conflictiva"
  )

# Identificar seeds (organizaciones políticas de pobladores, tipo 3)
seeds <- which(attr$tipo == "3")

# Verificar que los datos se cargaron correctamente
print("Dimensiones de las redes:")
print(dim(conf_net))
print("Número de organizaciones semilla:")
print(length(seeds))
print("Nombres de las organizaciones semilla:")
print(attr$Nombre[seeds])

# Verificar organizaciones conflictivas iniciales
print("Número de organizaciones inicialmente conflictivas:")
print(sum(attr$conflictivo))





improved_diffusion_v4 <- function(networks, attributes, seeds,
                                  thresholds = c(0.08, 0.12, 0.15, 0.2)) {
  
  # Definir pesos según teoría de movimientos sociales
  # Damos más peso a recursos y valores compartidos que reflejan mejor 
  # la coordinación del movimiento social (Diani, 2015)
  weights <- c(0.25, 0.35, 0.15, 0.25) # confianza, valores, cooperación, recursos
  
  n_nodes <- nrow(networks[[1]])
  composite_net <- matrix(0, n_nodes, n_nodes)
  
  for(i in 1:length(networks)) {
    composite_net <- composite_net + networks[[i]] * weights[i]
  }
  
  # Identificar roles organizativos clave basados en Diani y Tilly
  g <- graph_from_adjacency_matrix(composite_net, weighted=TRUE) 
  
  # Identificar brokers usando múltiples medidas
  between <- betweenness(g, normalized=TRUE)
  flow_bet <- edge_betweenness(g)
  deg <- degree(g, normalized=TRUE)
  eigen <- eigen_centrality(g)$vector
  
  # Calcular score de influencia organizativa
  org_influence <- scale(between) + scale(deg) + scale(eigen)
  
  # Identificar organizaciones clave según su posición estructural
  key_orgs <- list(
    brokers = which(org_influence > quantile(org_influence, 0.8)),
    movement_cores = which(deg > quantile(deg, 0.8)),
    bridges = which(between > quantile(between, 0.8))
  )
  
  # Estados iniciales 
  states <- rep(0, n_nodes)
  states[seeds] <- 1
  states[which(attributes$conflictivo)] <- 1
  
  history <- matrix(0, nrow=50, ncol=n_nodes)
  history[1,] <- states
  
  # Proceso de difusión mejorado
  for(iter in 2:50) {
    old_states <- states
    
    for(i in which(states == 0)) {
      
      # Calcular influencia por tipo de vínculo
      influence <- 0
      total_weight <- 0
      
      for(n in 1:length(networks)) {
        neighbors <- which(networks[[n]][i,] > 0)
        if(length(neighbors) > 0) {
          # Influencia ponderada por tipo de vínculo
          net_influence <- sum(states[neighbors]) / length(neighbors)
          
          # Bonus por tipo de vínculo según teoría
          if(n == 2) { # Valores compartidos
            net_influence <- net_influence * 1.2
          }
          if(n == 4) { # Recursos compartidos  
            net_influence <- net_influence * 1.1
          }
          
          influence <- influence + net_influence * weights[n]
          total_weight <- total_weight + weights[n]
        }
      }
      
      if(total_weight > 0) {
        influence <- influence / total_weight
        
        # Mecanismos de activación de identidades políticas
        
        # 1. Efecto de brokers (Tilly)
        if(any(neighbors %in% key_orgs$brokers)) {
          influence <- influence * 1.3
        }
        
        # 2. Efecto de organizaciones núcleo del movimiento (Diani)
        if(any(neighbors %in% key_orgs$movement_cores)) {
          influence <- influence * 1.2  
        }
        
        # 3. Bonus por multiplexidad (redes múltiples indican mayor compromiso)
        multiplex_bonus <- 0
        for(j in neighbors) {
          ties <- sum(sapply(networks, function(x) x[i,j] > 0))
          if(ties > 1) multiplex_bonus <- multiplex_bonus + 0.15 * (ties-1)
        }
        influence <- influence + multiplex_bonus
        
        # 4. Efecto del bloque estructural (posición en la red)
        block <- attributes$bloque_concor[i]
        if(!is.na(block)) {
          threshold <- thresholds[block]
          # Ajustar umbral según tipo organizativo
          if(attributes$tipo[i] %in% c(2,3,4)) { # OPP, CV y OPB
            threshold <- threshold * 0.8 # Más susceptibles a activación
          }
          if(influence >= threshold) {
            states[i] <- 1
          }
        }
      }
    }
    
    history[iter,] <- states
    if(all(old_states == states)) break
  }
  
  return(list(
    final_states = states,
    history = history[1:iter,],
    n_iterations = iter,
    key_orgs = key_orgs,
    composite_net = composite_net,
    influence_scores = org_influence
  ))
}



analyze_diffusion_results <- function(results, attributes) {
  
  # Análisis básico de difusión
  diffusion_df <- data.frame(
    nombre = attributes$Nombre,
    tipo = attributes$tipo, 
    bloque = attributes$bloque_concor,
    is_broker = 1:nrow(attributes) %in% results$key_orgs$brokers,
    is_core = 1:nrow(attributes) %in% results$key_orgs$movement_cores,
    is_bridge = 1:nrow(attributes) %in% results$key_orgs$bridges,
    influence_score = results$influence_scores,
    inicial = 1:nrow(attributes) %in% seeds | attributes$conflictivo,
    final = as.logical(results$final_states)
  )
  
  # Análisis por tipo organizativo
  type_results <- diffusion_df %>%
    group_by(tipo) %>%
    summarise(
      n = n(),
      inicial = sum(inicial),
      final = sum(final), 
      tasa_adopcion = (sum(final) - sum(inicial))/(n - sum(inicial)),
      n_brokers = sum(is_broker),
      n_cores = sum(is_core),
      influence_mean = mean(influence_score)
    )
  
  # Análisis por bloque estructural  
  block_results <- diffusion_df %>%
    group_by(bloque) %>%
    summarise(
      n = n(),
      inicial = sum(inicial),
      final = sum(final),
      tasa_adopcion = (sum(final) - sum(inicial))/(n - sum(inicial)),
      n_brokers = sum(is_broker), 
      n_cores = sum(is_core),
      influence_mean = mean(influence_score)
    )
  
  # Identificar organizaciones clave
  key_orgs <- diffusion_df %>%
    filter(is_broker | is_core | is_bridge) %>%
    arrange(desc(influence_score)) %>%
    select(nombre, tipo, bloque, influence_score,
           is_broker, is_core, is_bridge)
  
  return(list(
    diffusion_df = diffusion_df,
    type_results = type_results, 
    block_results = block_results,
    key_orgs = key_orgs
  ))
}





# Ejecutar el modelo mejorado
results <- improved_diffusion_v4(
  networks = list(conf_net, val_net, coop_net, rec_net),
  attributes = attr,
  seeds = seeds
)

# Analizar resultados
analysis <- analyze_diffusion_results(results, attr)

# Ver resultados principales
print("Resultados por tipo de organización:")
print(analysis$type_results)

print("\nResultados por bloque estructural:")
print(analysis$block_results)

print("\nOrganizaciones clave en el proceso:")
print(analysis$key_orgs)

# Visualizar la evolución de la difusión
plot_diffusion <- function(history) {
  adoption_over_time <- colSums(history) / ncol(history)
  plot(adoption_over_time, type = "l", 
       xlab = "Iteración", 
       ylab = "Proporción de organizaciones activadas",
       main = "Evolución de la difusión")
}

plot_diffusion(results$history)


# Análisis de la difusión por tipo de organización
detailed_type_analysis <- analysis$diffusion_df %>%
  group_by(tipo) %>%
  summarise(
    total = n(),
    inicial_conflictivas = sum(inicial),
    final_conflictivas = sum(final),
    tasa_adopcion = (sum(final) - sum(inicial))/(n() - sum(inicial)),
    prop_brokers = mean(is_broker),
    prop_core = mean(is_core),
    influence_promedio = mean(influence_score)
  ) %>%
  arrange(desc(tasa_adopcion))

print("Análisis detallado por tipo de organización:")
print(detailed_type_analysis)

# Identificar las organizaciones más influyentes
top_influential <- analysis$diffusion_df %>%
  arrange(desc(influence_score)) %>%
  select(nombre, tipo, influence_score, is_broker, is_core, is_bridge) %>%
  head(10)

print("\nLas 10 organizaciones más influyentes:")
print(top_influential)

# Análisis de la multiplexidad de vínculos
calc_multiplexity <- function(networks, node) {
  ties <- sapply(networks, function(x) which(x[node,] > 0))
  multiplex_ties <- Reduce(intersect, ties)
  return(length(multiplex_ties))
}

multiplexity_scores <- sapply(1:nrow(attr), function(i) {
  calc_multiplexity(list(conf_net, val_net, coop_net, rec_net), i)
})

print("\nOrganizaciones con mayor multiplexidad de vínculos:")
data.frame(
  nombre = attr$Nombre,
  multiplexidad = multiplexity_scores
) %>%
  arrange(desc(multiplexidad)) %>%
  head(10) %>%
  print()









