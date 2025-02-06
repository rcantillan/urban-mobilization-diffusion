---
title: "Análisis de Difusión de Identidades Políticas en Redes de Organizaciones"
subtitle: "Un estudio del conflicto por el Plan Regulador de Peñalolén 2011"
author: "Roberto Cantillan"
draft: true
---

## 1. Introducción

### 1.1 Contexto del Estudio
El conflicto por el Plan Regulador Comunal de Peñalolén en 2011 representa un caso paradigmático de movilización social en el Chile neoliberal. Este episodio permite analizar cómo emergen identidades políticas a través de la activación de redes organizativas preexistentes.

### 1.2 Marco Teórico
La investigación se fundamenta en dos tradiciones teóricas complementarias:

1. **Teoría de Movimientos Sociales** (Diani, 2015):
   - Los movimientos sociales son modos de coordinación de acción colectiva
   - Se caracterizan por densas redes informales entre organizaciones
   - Comparten una identidad colectiva que trasciende eventos específicos

2. **Teoría de la Contienda Política** (Tilly, McAdam, Tarrow):
   - Énfasis en mecanismos y procesos causales
   - Importancia de la activación de fronteras identitarias
   - Rol de brokers en la difusión de marcos interpretativos

### 1.3 Hipótesis

1. **H1:** Las organizaciones políticas de pobladores (OPP) actúan como núcleo del proceso de difusión debido a su alta centralidad y multiplexidad de vínculos.

2. **H2:** La difusión de identidades políticas es más efectiva a través de vínculos que combinan valores compartidos y recursos.

3. **H3:** Los bloques estructurales identificados mediante CONCOR corresponden a dominios de acción colectiva diferenciados pero interconectados.

## 2. Datos y Métodos

### 2.1 Datos
```{r}
#| label: setup
#| message: false

# Cargar librerías
library(tidyverse)
library(igraph)
library(ggplot2)
library(gridExtra)
library(sna)
library(network)

# Cargar las matrices de redes desde los archivos CSV
conf_net <- as.matrix(read.csv("datos/ConfianzaNet.csv", header=FALSE))
coop_net <- as.matrix(read.csv("datos/CoopNet.csv", header=FALSE))
rec_net <-  as.matrix(read.csv("datos/RecursosNet.csv", header=FALSE))
val_net <-  as.matrix(read.csv("datos/ValoresNet.csv", header=FALSE))

# Cargar atributos
attr <- read.csv("datos/Atributos _org2011.csv") %>%
  mutate(
    bloque_concor = X1PosciónCONCOR,
    tipo = Tipo,
    conflictivo = trimws(Orientación) == "Conflictiva"
  )

# Identificar seeds
seeds <- which(attr$tipo == "3")
```


```{r}
#| label: network-description

# Función para calcular estadísticas descriptivas de redes
analyze_networks <- function(networks, names) {
  map2_dfr(networks, names, ~{
    tibble(
      red = .y,
      densidad = sum(.x > 0)/(nrow(.x)^2),
      grado_medio = mean(colSums(.x > 0)),
      reciprocidad = reciprocity(graph_from_adjacency_matrix(.x))
    )
  })
}

network_stats <- analyze_networks(
  list(conf_net, val_net, coop_net, rec_net),
  c("Confianza", "Valores", "Cooperación", "Recursos")
)

print(network_stats)
```


```{r}
improved_diffusion_v5 <- function(networks, attributes, seeds,
                                max_iterations = 100,
                                convergence_threshold = 0.01,
                                thresholds = c(0.08, 0.12, 0.15, 0.2)) {
  
  # Pesos de las redes según teoría de movimientos sociales (Diani)
  # Valores y confianza son más importantes para identidad colectiva
  # Recursos y cooperación para capacidad de movilización
  weights <- c(0.25, 0.35, 0.15, 0.25) # confianza, valores, cooperación, recursos
  
  n_nodes <- nrow(networks[[1]])
  composite_net <- matrix(0, n_nodes, n_nodes)
  
  # Crear red compuesta ponderada
  for(i in 1:length(networks)) {
    composite_net <- composite_net + networks[[i]] * weights[i]
  }
  
  # Crear grafo para medidas de centralidad
  g <- graph_from_adjacency_matrix(composite_net, weighted=TRUE) 
  
  # Calcular múltiples medidas de centralidad
  between <- igraph::betweenness(g, normalized=TRUE)
  flow_bet <- edge_betweenness(g)
  deg <- igraph::degree(g, normalized=TRUE)
  eigen <- eigen_centrality(g)$vector
  close <- igraph::closeness(g, normalized=TRUE)
  
  # Score compuesto de influencia organizativa
  # Incorpora múltiples dimensiones de centralidad
  org_influence <- scale(between) + 
                  scale(deg) + 
                  scale(eigen) + 
                  scale(close)
  
  # Identificar roles organizativos clave (Diani & Tilly)
  key_orgs <- list(
    # Brokers: alto between y flow
    brokers = which(between > quantile(between, 0.8)),
    
    # Cores: alto degree y eigenvector
    movement_cores = which(deg > quantile(deg, 0.8)),
    
    # Bridges: alto betweenness pero bajo degree
    bridges = which(between > quantile(between, 0.8) & 
                   deg < median(deg))
  )
  
  # Calcular multiplexidad de vínculos
  calc_multiplexity <- function(node) {
    ties <- sapply(networks, function(x) which(x[node,] > 0))
    multiplex_ties <- Reduce(intersect, ties)
    return(length(multiplex_ties))
  }
  multiplexity <- sapply(1:n_nodes, calc_multiplexity)
  
  # Estados iniciales
  states <- rep(0, n_nodes)
  states[seeds] <- 1  # Seeds son OPP
  states[which(attributes$conflictivo)] <- 1  # Ya conflictivos
  
  # Matriz para guardar historia del proceso
  history <- matrix(0, nrow=max_iterations, ncol=n_nodes)
  history[1,] <- states
  
  # Variables para convergencia
  recent_changes <- numeric(5)
  
  # Proceso de difusión
  for(iter in 2:max_iterations) {
    old_states <- states
    
    # Para cada nodo no activado
    for(i in which(states == 0)) {
      
      # Calcular influencia por tipo de vínculo
      influence <- 0
      total_weight <- 0
      
      for(n in 1:length(networks)) {
        neighbors <- which(networks[[n]][i,] > 0)
        if(length(neighbors) > 0) {
          # Influencia base: proporción de vecinos activados
          net_influence <- sum(states[neighbors]) / length(neighbors)
          
          # Bonus por tipo de vínculo
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
        
        # Mecanismos de activación basados en teoría
        
        # 1. Efecto de brokers (Tilly)
        if(any(neighbors %in% key_orgs$brokers)) {
          influence <- influence * 1.3
        }
        
        # 2. Efecto de organizaciones núcleo (Diani)
        if(any(neighbors %in% key_orgs$movement_cores)) {
          influence <- influence * 1.2
        }
        
        # 3. Efecto de puentes estructurales
        if(any(neighbors %in% key_orgs$bridges)) {
          influence <- influence * 1.1
        }
        
        # 4. Bonus por multiplexidad 
        multiplex_bonus <- 0.15 * multiplexity[i]
        influence <- influence + multiplex_bonus
        
        # Ajuste de umbral según tipo organizativo
        block <- attributes$bloque_concor[i]
        if(!is.na(block)) {
          threshold <- thresholds[block]
          
          # Tipos más susceptibles a activación
          if(attributes$tipo[i] %in% c(2,3,4)) { # OPP, CV y OPB
            threshold <- threshold * 0.8
          }
          
          # Verificar activación
          if(influence >= threshold) {
            states[i] <- 1
          }
        }
      }
    }
    
    # Guardar historia
    history[iter,] <- states
    
    # Verificar convergencia
    if(all(old_states == states)) break
    
    # Calcular tasa de cambio reciente
    if(iter > 5) {
      recent_changes <- c(recent_changes[-1], 
                         mean(states) - mean(old_states))
      if(abs(mean(recent_changes)) < convergence_threshold) break
    }
  }
  
  # Recortar matriz de historia
  history <- history[1:iter,]
  
  return(list(
    final_states = states,
    history = history,
    n_iterations = iter,
    key_orgs = key_orgs,
    composite_net = composite_net,
    influence_scores = org_influence,
    multiplexity = multiplexity,
    converged = iter < max_iterations,
    final_change_rate = mean(recent_changes)
  ))
}
```





```{r}
plot_diffusion_improved <- function(history) {
  # Asegurarnos que usamos las dimensiones correctas
  # history es una matriz donde las filas son iteraciones y columnas son nodos
  prop_activated <- rowMeans(history)  # Promedio por fila en lugar de colSums
  
  # Crear dataframe para ggplot
  df <- data.frame(
    iteracion = seq_len(length(prop_activated)),
    proporcion = prop_activated
  )
  
  # Crear visualización
  ggplot(df, aes(x = iteracion, y = proporcion)) +
    geom_line(color = "#2C3E50", size = 1) +
    geom_point(color = "#E74C3C", size = 2) +
    theme_minimal() +
    labs(
      title = "Evolución del Proceso de Difusión",
      subtitle = paste("Proceso convergió en", nrow(history), "iteraciones"),
      x = "Iteración",
      y = "Proporción de organizaciones activadas"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12, color = "grey40"),
      axis.title = element_text(size = 11),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 1)
    ) +
    geom_hline(
      yintercept = max(prop_activated), 
      linetype = "dashed", 
      color = "grey50"
    ) +
    annotate(
      "text",
      x = max(df$iteracion)/2,
      y = max(prop_activated) + 0.05,
      label = sprintf("Máximo: %.1f%%", max(prop_activated)*100),
      color = "grey40"
    )
}

plot_diffusion_detailed <- function(history, attributes) {
  # Calcular proporciones por tipo de organización
  prop_by_type <- apply(history, 1, function(x) {
    tapply(x, attributes$tipo, mean)
  }) %>% t() %>% as.data.frame()
  
  # Añadir iteración
  prop_by_type$iteracion <- 1:nrow(prop_by_type)
  
  # Convertir a formato largo
  df_long <- prop_by_type %>%
    pivot_longer(
      cols = -iteracion,
      names_to = "tipo",
      values_to = "proporcion"
    )
  
  # Crear visualización
  ggplot(df_long, aes(x = iteracion, y = proporcion, color = tipo)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(
      title = "Evolución de la Difusión por Tipo de Organización",
      subtitle = paste("Proceso convergió en", nrow(history), "iteraciones"),
      x = "Iteración",
      y = "Proporción activada",
      color = "Tipo de organización"
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_brewer(palette = "Set2") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

```



```{r}
#| label: diffusion-patterns

# Ejecutar modelo
results <- improved_diffusion_v5(
  networks = list(conf_net, val_net, coop_net, rec_net),
  attributes = attr,
  seeds = seeds
)

# Visualizar evolución
plot_diffusion_improved(results$history)
```



```{r}
# Visualizar ambos plots juntos
library(gridExtra)
grid.arrange(
  plot_diffusion_improved(results$history),
  plot_diffusion_detailed(results$history, attr),
  ncol = 2
)

```


```{r}
#| label: block-analysis
#| eval: false

# Visualizar resultados por bloque
ggplot(analysis$block_results, 
       aes(x = as.factor(bloque), 
           y = tasa_adopcion)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Tasa de Adopción por Bloque Estructural")

# Análisis de composición de bloques
block_composition <- analysis$diffusion_df %>%
  group_by(bloque) %>%
  count(tipo) %>%
  spread(tipo, n, fill = 0)

print(block_composition)
```


```{r}
#| label: key-orgs
#| eval: false

# Identificar y visualizar organizaciones clave
plot_network_improved(g_composite, attr)

# Análisis de multiplexidad
print(analysis$key_orgs)
```







