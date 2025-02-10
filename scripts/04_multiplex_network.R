##############################################################################
# EJEMPLO COMPLETO EN R: 3D MULTIPLEX con INTRA-layer en NEGRO y
# INTER-layer en ROSA, mostrando nombres de cada layer
# y una nota/footnote en inglés sobre la definición de "inter-layer link".
##############################################################################

library(igraph)
library(plotly)

##############################################################################
# 1) CARGAR MATRICES (3 capas de ejemplo)
##############################################################################
coop_matrix  <- as.matrix(read.csv("data/CoopNet.csv",       header=FALSE))
trust_matrix <- as.matrix(read.csv("data/ConfianzaNet.csv",  header=FALSE))
reso_matrix  <- as.matrix(read.csv("data/RecursosNet.csv",   header=FALSE))

n <- nrow(coop_matrix)

g_coop  <- graph_from_adjacency_matrix(coop_matrix,  mode="directed")
g_trust <- graph_from_adjacency_matrix(trust_matrix, mode="directed")
g_reso  <- graph_from_adjacency_matrix(reso_matrix,  mode="directed")

node_names <- paste0("Node", 1:n)
V(g_coop)$name  <- node_names
V(g_trust)$name <- node_names
V(g_reso)$name  <- node_names

##############################################################################
# 2) EXTRAER ARISTAS INTRALAYER
##############################################################################
extract_intralayer_edges <- function(adj_matrix, layer_id) {
  df <- data.frame()
  N  <- nrow(adj_matrix)
  for(i in 1:N) {
    for(j in 1:N) {
      if(adj_matrix[i,j] != 0) {
        df <- rbind(df, data.frame(from=i, to=j, layer=layer_id))
      }
    }
  }
  return(df)
}

edges_coop  <- extract_intralayer_edges(coop_matrix,  1)
edges_trust <- extract_intralayer_edges(trust_matrix, 2)
edges_reso  <- extract_intralayer_edges(reso_matrix,  3)
edges_intralayer <- rbind(edges_coop, edges_trust, edges_reso)

##############################################################################
# 3) CREAR ARISTAS INTERLAYER (p.ej., cuando (i->j) aparece en más de 1 capa)
##############################################################################
create_interlayer_edges <- function(adj_mat_A, layerA, adj_mat_B, layerB) {
  df <- data.frame()
  N  <- nrow(adj_mat_A)
  for(i in 1:N) {
    for(j in 1:N) {
      if(i != j) {
        if(adj_mat_A[i,j] != 0 && adj_mat_B[i,j] != 0) {
          df <- rbind(df, data.frame(
            from_id    = i,
            from_layer = layerA,
            to_id      = j,
            to_layer   = layerB
          ))
        }
      }
    }
  }
  return(df)
}

edges_inter_coop_trust <- create_interlayer_edges(coop_matrix,  1,
                                                  trust_matrix, 2)
edges_inter_coop_reso  <- create_interlayer_edges(coop_matrix,  1,
                                                  reso_matrix,  3)
edges_inter_trust_reso <- create_interlayer_edges(trust_matrix, 2,
                                                  reso_matrix,  3)

edges_interlayer <- rbind(edges_inter_coop_trust,
                          edges_inter_coop_reso,
                          edges_inter_trust_reso)

##############################################################################
# 4) DEFINIR LAS POSICIONES (x,y) Y REPLICAR PARA z=1,2,3
##############################################################################
layout_2d <- layout_with_fr(g_coop)
layout_2d <- as.data.frame(layout_2d)
colnames(layout_2d) <- c("x","y")

df_nodes_layer1 <- data.frame(node_id=1:n, node_name=node_names, layer=1,
                              x=layout_2d$x, y=layout_2d$y, z=1)
df_nodes_layer2 <- data.frame(node_id=1:n, node_name=node_names, layer=2,
                              x=layout_2d$x, y=layout_2d$y, z=2)
df_nodes_layer3 <- data.frame(node_id=1:n, node_name=node_names, layer=3,
                              x=layout_2d$x, y=layout_2d$y, z=3)

df_nodes_all <- rbind(df_nodes_layer1, df_nodes_layer2, df_nodes_layer3)

##############################################################################
# 5) PLANOS DE CAPA (colores semitransparentes)
##############################################################################
xrange <- range(df_nodes_all$x)
yrange <- range(df_nodes_all$y)

add_layer_plane <- function(plot_obj, z_val, plane_color="gray", alpha=0.2) {
  x_corners <- c(xrange[1], xrange[2], xrange[2], xrange[1])
  y_corners <- c(yrange[1], yrange[1], yrange[2], yrange[2])
  z_corners <- rep(z_val, 4)
  
  i_ <- c(0,0)
  j_ <- c(1,2)
  k_ <- c(2,3)
  
  plot_obj %>%
    add_trace(
      type="mesh3d",
      x = x_corners,
      y = y_corners,
      z = z_corners,
      i = i_, j = j_, k = k_,
      opacity = alpha,
      facecolor = rep(plane_color, length(i_)),
      hoverinfo = "none",
      showlegend=FALSE
    )
}

##############################################################################
# 6) OPCIONAL: ELIMINAR NODOS QUE NO TIENEN VINCULOS
##############################################################################
# Filtra nodos sin lazos ni intra-layer ni inter-layer
active_nodes <- unique(c(
  edges_intralayer$from, 
  edges_intralayer$to,
  edges_interlayer$from_id, 
  edges_interlayer$to_id
))
df_nodes_all <- subset(df_nodes_all, node_id %in% active_nodes)
edges_intralayer <- subset(
  edges_intralayer,
  from %in% active_nodes & to %in% active_nodes
)
edges_interlayer <- subset(
  edges_interlayer,
  from_id %in% active_nodes & to_id %in% active_nodes
)

##############################################################################
# 7) GRAFICAR EN 3D: Aristas intra (negro) e inter (rosa) + Nombre de ejes
##############################################################################
p <- plot_ly()

# (A) Dibujamos un plano para cada capa con color distinto
p <- add_layer_plane(p, z_val=1, plane_color="red",   alpha=0.2)
p <- add_layer_plane(p, z_val=2, plane_color="green", alpha=0.2)
p <- add_layer_plane(p, z_val=3, plane_color="blue",  alpha=0.2)

# (B) Nodos
p <- p %>%
  add_markers(
    data = df_nodes_all,
    x=~x, y=~y, z=~z,
    marker=list(size=6, color=~factor(layer), colorscale="Viridis"),
    text=~paste("Node:", node_name, "<br>Layer:", layer),
    hoverinfo="text",
    name="Nodos"
  )

# (C) Aristas INTRA-layer (negro)
for(k in 1:nrow(edges_intralayer)) {
  lay <- edges_intralayer$layer[k]
  fn  <- edges_intralayer$from[k]
  tn  <- edges_intralayer$to[k]
  
  coords_from <- subset(df_nodes_all, node_id==fn & layer==lay)
  coords_to   <- subset(df_nodes_all, node_id==tn & layer==lay)
  
  # >>> AQUI: COLOR, ANCHO Y TRANSPARENCIA DE ARISTAS INTRA (negro)
  edge_color <- "black"
  edge_width <- 1.8
  edge_alpha <- 0.8  # no transparency
  
  p <- p %>%
    add_trace(
      x = c(coords_from$x, coords_to$x, NA),
      y = c(coords_from$y, coords_to$y, NA),
      z = c(coords_from$z, coords_to$z, NA),
      mode="lines",
      line=list(color=edge_color, width=edge_width),
      opacity=edge_alpha,
      hoverinfo="none",
      showlegend=FALSE
    )
}

# (D) Aristas INTER-layer (rosa)
for(k in 1:nrow(edges_interlayer)) {
  f_id    <- edges_interlayer$from_id[k]
  f_layer <- edges_interlayer$from_layer[k]
  t_id    <- edges_interlayer$to_id[k]
  t_layer <- edges_interlayer$to_layer[k]
  
  coords_from <- subset(df_nodes_all, node_id==f_id & layer==f_layer)
  coords_to   <- subset(df_nodes_all, node_id==t_id & layer==t_layer)
  
  # >>> AQUI: COLOR, ANCHO Y TRANSPARENCIA DE ARISTAS INTER (rosa)
  edge_color <- "red"
  edge_width <- 1
  edge_alpha <- 0.5  # no transparency
  
  p <- p %>%
    add_trace(
      x = c(coords_from$x, coords_to$x, NA),
      y = c(coords_from$y, coords_to$y, NA),
      z = c(coords_from$z, coords_to$z, NA),
      mode="lines",
      line=list(color=edge_color, width=edge_width, dash="solid"),
      opacity=edge_alpha,
      hoverinfo="none",
      showlegend=FALSE
    )
}



p <- p %>% layout(
  title = "Multiplex Network (Intra=Black, Inter=Pink)",
  scene = list(
    xaxis = list(title=""),
    yaxis = list(title=""),
    zaxis = list(
      title="", 
      tickmode="array",
      tickvals=c(1,2,3),
      ticktext=c("Cooperation", "Trust", "Resources")
    )
  ),
  annotations = list(
    x = 1,  # posición en el plano
    y = 0,  # posición en el plano
    text = "",
    xref = "paper",
    yref = "paper",
    showarrow = FALSE,
    xanchor = "right",
    yanchor = "auto",
    font = list(size=10, color="gray40")
  )
)

# Mostramos
p

##############################################################################
# Created/Modified: ["data/CoopNet.csv", "data/ConfianzaNet.csv","data/RecursosNet.csv"]
##############################################################################
library(igraph)
library(plotly)

##############################################################################
# 1) CARGAR MATRICES (5 capas)
# Cada matriz representa una capa diferente de la red multiplex
##############################################################################
coop_matrix    <- as.matrix(read.csv("data/CoopNet.csv",       header=FALSE))
trust_matrix   <- as.matrix(read.csv("data/ConfianzaNet.csv",  header=FALSE))
reso_matrix    <- as.matrix(read.csv("data/RecursosNet.csv",   header=FALSE))
values_matrix  <- as.matrix(read.csv("data/ValoresNet.csv",    header=FALSE))
kinship_matrix <- as.matrix(read.csv("data/ParentescoNet.csv", header=FALSE))

n <- nrow(coop_matrix)

# Convertir matrices a objetos igraph para facilitar el manejo
g_coop    <- graph_from_adjacency_matrix(coop_matrix,    mode="directed")
g_trust   <- graph_from_adjacency_matrix(trust_matrix,   mode="directed")
g_reso    <- graph_from_adjacency_matrix(reso_matrix,    mode="directed")
g_values  <- graph_from_adjacency_matrix(values_matrix,  mode="directed")
g_kinship <- graph_from_adjacency_matrix(kinship_matrix, mode="directed")

# Asignar nombres a los nodos para mejor identificación
node_names <- paste0("Node", 1:n)
V(g_coop)$name    <- node_names
V(g_trust)$name   <- node_names
V(g_reso)$name    <- node_names
V(g_values)$name  <- node_names
V(g_kinship)$name <- node_names

##############################################################################
# 2) EXTRAER ARISTAS INTRALAYER (conexiones dentro de cada capa)
##############################################################################
extract_intralayer_edges <- function(adj_matrix, layer_id) {
  df <- data.frame()
  N  <- nrow(adj_matrix)
  for(i in 1:N) {
    for(j in 1:N) {
      if(adj_matrix[i,j] != 0) {
        df <- rbind(df, data.frame(from=i, to=j, layer=layer_id))
      }
    }
  }
  return(df)
}

edges_coop    <- extract_intralayer_edges(coop_matrix,    1)
edges_trust   <- extract_intralayer_edges(trust_matrix,   2)
edges_reso    <- extract_intralayer_edges(reso_matrix,    3)
edges_values  <- extract_intralayer_edges(values_matrix,  4)
edges_kinship <- extract_intralayer_edges(kinship_matrix, 5)

edges_intralayer <- rbind(edges_coop, edges_trust, edges_reso, 
                          edges_values, edges_kinship)

##############################################################################
# 3) CREAR ARISTAS INTERLAYER (conexiones entre capas)
##############################################################################
create_interlayer_edges <- function(adj_mat_A, layerA, adj_mat_B, layerB) {
  df <- data.frame()
  N  <- nrow(adj_mat_A)
  for(i in 1:N) {
    for(j in 1:N) {
      if(i != j) {
        if(adj_mat_A[i,j] != 0 && adj_mat_B[i,j] != 0) {
          df <- rbind(df, data.frame(
            from_id    = i,
            from_layer = layerA,
            to_id      = j,
            to_layer   = layerB
          ))
        }
      }
    }
  }
  return(df)
}

# Crear todas las relaciones posibles entre capas
edges_inter_coop_trust <- create_interlayer_edges(coop_matrix,    1, trust_matrix,   2)
edges_inter_coop_reso  <- create_interlayer_edges(coop_matrix,    1, reso_matrix,    3)
edges_inter_coop_val   <- create_interlayer_edges(coop_matrix,    1, values_matrix,  4)
edges_inter_coop_kin   <- create_interlayer_edges(coop_matrix,    1, kinship_matrix, 5)
edges_inter_trust_reso <- create_interlayer_edges(trust_matrix,   2, reso_matrix,    3)
edges_inter_trust_val  <- create_interlayer_edges(trust_matrix,   2, values_matrix,  4)
edges_inter_trust_kin  <- create_interlayer_edges(trust_matrix,   2, kinship_matrix, 5)
edges_inter_reso_val   <- create_interlayer_edges(reso_matrix,    3, values_matrix,  4)
edges_inter_reso_kin   <- create_interlayer_edges(reso_matrix,    3, kinship_matrix, 5)
edges_inter_val_kin    <- create_interlayer_edges(values_matrix,  4, kinship_matrix, 5)

edges_interlayer <- rbind(
  edges_inter_coop_trust, edges_inter_coop_reso,  edges_inter_coop_val,
  edges_inter_coop_kin,   edges_inter_trust_reso, edges_inter_trust_val,
  edges_inter_trust_kin,  edges_inter_reso_val,   edges_inter_reso_kin,
  edges_inter_val_kin
)

##############################################################################
# 4) DEFINIR POSICIONES (x,y) Y REPLICAR PARA z CON MAYOR SEPARACIÓN
##############################################################################
layout_2d <- layout_with_fr(g_coop)
layout_2d <- as.data.frame(layout_2d)
colnames(layout_2d) <- c("x","y")

# Aumentamos la separación entre capas usando z=1,3,5,7,9
df_nodes_layer1 <- data.frame(node_id=1:n, node_name=node_names, layer=1,
                              x=layout_2d$x, y=layout_2d$y, z=2)
df_nodes_layer2 <- data.frame(node_id=1:n, node_name=node_names, layer=2,
                              x=layout_2d$x, y=layout_2d$y, z=8)
df_nodes_layer3 <- data.frame(node_id=1:n, node_name=node_names, layer=3,
                              x=layout_2d$x, y=layout_2d$y, z=14)
df_nodes_layer4 <- data.frame(node_id=1:n, node_name=node_names, layer=4,
                              x=layout_2d$x, y=layout_2d$y, z=20)
df_nodes_layer5 <- data.frame(node_id=1:n, node_name=node_names, layer=5,
                              x=layout_2d$x, y=layout_2d$y, z=26)

df_nodes_all <- rbind(df_nodes_layer1, df_nodes_layer2, df_nodes_layer3,
                      df_nodes_layer4, df_nodes_layer5)

##############################################################################
# 5) PLANOS DE CAPA (colores semitransparentes)
##############################################################################
xrange <- range(df_nodes_all$x)
yrange <- range(df_nodes_all$y)

add_layer_plane <- function(plot_obj, z_val, plane_color="gray", alpha=0.2) {
  x_corners <- c(xrange[1], xrange[2], xrange[2], xrange[1])
  y_corners <- c(yrange[1], yrange[1], yrange[2], yrange[2])
  z_corners <- rep(z_val, 4)
  
  i_ <- c(0,0)
  j_ <- c(1,2)
  k_ <- c(2,3)
  
  plot_obj %>%
    add_trace(
      type="mesh3d",
      x = x_corners,
      y = y_corners,
      z = z_corners,
      i = i_, j = j_, k = k_,
      opacity = alpha,
      facecolor = rep(plane_color, length(i_)),
      hoverinfo = "none",
      showlegend=FALSE
    )
}

##############################################################################
# 6) GRAFICAR EN 3D
##############################################################################
p <- plot_ly()

# Planos para cada capa con sus nuevas posiciones z
p <- add_layer_plane(p, z_val=2,  plane_color="red",    alpha=0.2)
p <- add_layer_plane(p, z_val=8,  plane_color="green",  alpha=0.2)
p <- add_layer_plane(p, z_val=14, plane_color="orange",   alpha=0.2)
p <- add_layer_plane(p, z_val=20, plane_color="purple", alpha=0.2)
p <- add_layer_plane(p, z_val=26, plane_color="blue", alpha=0.2)

# Nodos
p <- p %>%
  add_markers(
    data = df_nodes_all,
    x=~x, y=~y, z=~z,
    marker=list(size=6, color=~factor(layer), colorscale="Viridis"),
    text=~paste("Node:", node_name, "<br>Layer:", layer),
    hoverinfo="text",
    name="Nodos"
  )

# Aristas INTRA-layer (negro)
for(k in 1:nrow(edges_intralayer)) {
  lay <- edges_intralayer$layer[k]
  fn  <- edges_intralayer$from[k]
  tn  <- edges_intralayer$to[k]
  
  coords_from <- subset(df_nodes_all, node_id==fn & layer==lay)
  coords_to   <- subset(df_nodes_all, node_id==tn & layer==lay)
  
  edge_color <- "black"
  edge_width <- 1.8
  edge_alpha <- 0.8
  
  p <- p %>%
    add_trace(
      x = c(coords_from$x, coords_to$x, NA),
      y = c(coords_from$y, coords_to$y, NA),
      z = c(coords_from$z, coords_to$z, NA),
      mode="lines",
      line=list(color=edge_color, width=edge_width),
      opacity=edge_alpha,
      hoverinfo="none",
      showlegend=FALSE
    )
}

# Aristas INTER-layer (rosa)
for(k in 1:nrow(edges_interlayer)) {
  f_id    <- edges_interlayer$from_id[k]
  f_layer <- edges_interlayer$from_layer[k]
  t_id    <- edges_interlayer$to_id[k]
  t_layer <- edges_interlayer$to_layer[k]
  
  coords_from <- subset(df_nodes_all, node_id==f_id & layer==f_layer)
  coords_to   <- subset(df_nodes_all, node_id==t_id & layer==t_layer)
  
  edge_color <- "red"
  edge_width <- 0.6
  edge_alpha <- 0.4
  
  p <- p %>%
    add_trace(
      x = c(coords_from$x, coords_to$x, NA),
      y = c(coords_from$y, coords_to$y, NA),
      z = c(coords_from$z, coords_to$z, NA),
      mode="lines",
      line=list(color=edge_color, width=edge_width, dash="solid"),
      opacity=edge_alpha,
      hoverinfo="none",
      showlegend=FALSE
    )
}

# Configuración final con rotación de 180 grados
p <- p %>% layout(
  title = "",
  scene = list(
    xaxis = list(title=""),
    yaxis = list(title=""),
    zaxis = list(
      title="", 
      tickmode="array",
      tickvals=c(2,8,14,20,26),
      ticktext=c("Cooperation", "Trust", "Resources", "Values", "Kinship"),
      tickangle=0  # Este es el parámetro que gira el texto
      
    ),
    camera = list(
      eye = list(x=0.1, y=-2, z=0.1)  # Vista vertical
    )
  )
)

# Mostrar el plot
p

