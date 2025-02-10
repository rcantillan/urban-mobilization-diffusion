# Cargar librerías necesarias
library(ergm.multi)
library(network)
library(tidyverse)
library(Matrix)

# 1. Limpiar y preparar los atributos
# ===============================================
attributes <- read.csv("data/Atributos _org2011.csv", encoding = "UTF-8")

clean_attributes <- attributes %>%
  mutate(
    # Limpieza básica de columnas categóricas
    tipo = Tipo,  # ya está limpio
    posicion = X1PosciónCONCOR,  # ya está limpio
    tipo_org = str_trim(Tipo_de_organización),  # quitar espacios extra
    ubicacion = str_trim(Ubicación),  # quitar espacios extra
    orientacion = str_trim(Orientación),  # quitar espacios extra
    num_repertorios = NumRep,  # ya está limpio
    
    # Variables adicionales que podrían ser útiles para el análisis
    betweenness = nBetweenness,
    in_degree = NrmInDeg
  )

# 2. Cargar y preparar las matrices de redes
# ===============================================
mats_paths <- list(
  Cooperation = "data/CoopNet.csv",
  Trust = "data/ConfianzaNet.csv", 
  Resources = "data/RecursosNet.csv",
  Values = "data/ValoresNet.csv",
  Kinship = "data/ParentescoNet.csv"
)

# Función para leer matrices
read_network_matrix <- function(path) {
  if (!file.exists(path)) stop(paste("Archivo no encontrado:", path))
  as.matrix(read.csv(path, header = FALSE))
}

network_matrices <- lapply(mats_paths, read_network_matrix)

# 3. Crear objetos network con atributos
# ===============================================
network_objects <- lapply(network_matrices, function(mat) {
  net <- network(mat, directed = TRUE)
  
  # Agregar atributos de manera segura
  set.vertex.attribute(net, "tipo", clean_attributes$tipo)
  set.vertex.attribute(net, "posicion", clean_attributes$posicion)
  set.vertex.attribute(net, "tipo_org", clean_attributes$tipo_org)
  set.vertex.attribute(net, "ubicacion", clean_attributes$ubicacion)
  set.vertex.attribute(net, "orientacion", clean_attributes$orientacion)
  set.vertex.attribute(net, "num_repertorios", clean_attributes$num_repertorios)
  set.vertex.attribute(net, "betweenness", clean_attributes$betweenness)
  set.vertex.attribute(net, "in_degree", clean_attributes$in_degree)
  
  return(net)
})

# 5. Ajustar modelo ERGM multicapa
# ===============================================
# Este modelo está diseñado específicamente para analizar la movilización social
# 1. Primero mantengamos la preparación de datos igual hasta la creación del objeto Layer
my_layers <- Layer(
  cooperation = network_objects$Cooperation,
  trust = network_objects$Trust,
  resources = network_objects$Resources,
  values = network_objects$Values,
  kinship = network_objects$Kinship
)

# 2. Modelo ERGM multicapa simplificado
ergm_fit_simple <- ergm(my_layers ~ 
                          # Efectos básicos de densidad para cada capa
                          L(~edges, ~cooperation) +
                          L(~edges, ~trust) +
                          L(~edges, ~resources) +
                          L(~edges, ~values) +
                          L(~edges, ~kinship) +
                          
                          # Reciprocidad dentro de cada capa
                          L(~mutual, ~cooperation) +
                          L(~mutual, ~trust) +
                          L(~mutual, ~resources) +
                          L(~mutual, ~values) +
                          L(~mutual, ~kinship) +
                          
                          # Solo un efecto básico de homofilia por tipo de organización
                          L(~nodematch("tipo_org"), ~cooperation) +
                          
                          # Una interdependencia simple entre capas clave
                          L(~edges, ~cooperation&trust),
                        
                        control = control.ergm(
                          MCMLE.maxit = 25,
                          MCMC.burnin = 25000,
                          MCMC.samplesize = 2500,
                          parallel = 8,
                          parallel.type = "PSOCK",
                          MCMC.samplesize = 10000
                        )
)


# 6. Análisis de resultados
# ===============================================
summary(ergm_fit_simple)

# Guardar resultados
save(ergm_fit_simple, file = "ergm_multilayer_results.RData")



# Modelo ERGM expandido
ergm_fit_expanded <- ergm(my_layers ~ 
                            # Mantenemos los efectos básicos de densidad
                            L(~edges, ~cooperation) +
                            L(~edges, ~trust) +
                            L(~edges, ~resources) +
                            L(~edges, ~values) +
                            L(~edges, ~kinship) +
                            
                            # Mantenemos los efectos de reciprocidad
                            L(~mutual, ~cooperation) +
                            L(~mutual, ~trust) +
                            L(~mutual, ~resources) +
                            L(~mutual, ~values) +
                            L(~mutual, ~kinship) +
                            
                            # Efectos de homofilia - agregamos orientación
                            L(~nodematch("tipo_org"), ~cooperation) +
                            L(~nodematch("orientacion"), ~cooperation) + # Nuevo efecto
                            
                            # Interdependencias entre capas - expandimos las relaciones
                            L(~edges, ~cooperation&trust) +
                            L(~edges, ~cooperation&resources) + # Nueva interdependencia
                            L(~edges, ~trust&resources),        # Nueva interdependencia
                          
                          control = control.ergm(
                            MCMLE.maxit = 60,
                            MCMC.burnin = 25000,
                            parallel = 8,
                            parallel.type = "PSOCK",
                            MCMC.samplesize = 10000
                          )
)


summary(ergm_fit_expanded)

################################################################################
# Modelo ERGM complejo ---------------------------------------------------------
################################################################################

ergm_fit_complex <- ergm(my_layers ~ 
                           # Efectos básicos de densidad
                           L(~edges, ~cooperation) +
                           L(~edges, ~trust) +
                           L(~edges, ~resources) +
                           L(~edges, ~values) +
                           L(~edges, ~kinship) +
                           
                           # Efectos de reciprocidad
                           L(~mutual, ~cooperation) +
                           L(~mutual, ~trust) +
                           L(~mutual, ~resources) +
                           L(~mutual, ~values) +
                           L(~mutual, ~kinship) +
                           
                           # Efectos de homofilia
                           L(~nodematch("tipo_org"), ~cooperation) +
                           L(~nodematch("orientacion"), ~cooperation) +
                           
                           # Interdependencias entre capas centrales
                           L(~edges, ~cooperation&trust) +
                           L(~edges, ~cooperation&resources) +
                           L(~edges, ~trust&resources) +
                           
                           # Efectos de clustering
                           #L(~gwesp(0.5, fixed=TRUE), ~cooperation) +
                           #L(~gwesp(0.5, fixed=TRUE), ~trust) +
                           
                           # Interdependencias con valores y parentesco
                           L(~edges, ~cooperation&values) +
                           L(~edges, ~cooperation&kinship) +    # Nueva: cooperación basada en parentesco
                           L(~edges, ~trust&kinship),          # Nueva: confianza basada en parentesco
                         
                         control = control.ergm(
                           MCMC.burnin = 25000,
                           MCMC.samplesize = 5000,
                           parallel = 8,
                           parallel.type = "PSOCK"
                         )
)


# Resumen del modelo
summary(ergm_fit_complex)



# Modelo ERGM con complejidad moderada
ergm_fit_moderate <- ergm(my_layers ~ 
                            # Efectos básicos de densidad
                            L(~edges, ~cooperation) +
                            L(~edges, ~trust) +
                            L(~edges, ~resources) +
                            L(~edges, ~values) +
                            L(~edges, ~kinship) +
                            
                            # Efectos de reciprocidad
                            L(~mutual, ~cooperation) +
                            L(~mutual, ~trust) +
                            L(~mutual, ~resources) +
                            L(~mutual, ~values) +
                            L(~mutual, ~kinship) +
                            
                            # Efectos de homofilia fundamentales
                            L(~nodematch("tipo_org"), ~cooperation) +
                            L(~nodematch("orientacion"), ~cooperation) +
                            
                            # Interdependencias clave entre capas
                            L(~edges, ~cooperation&trust) +
                            L(~edges, ~cooperation&kinship) +
                            L(~edges, ~trust&kinship) +
                            
                            # Solo un efecto de clustering en cooperación
                            L(~gwesp(0.5, fixed=TRUE), ~cooperation),
                          
                          control = control.ergm(
                            MCMC.burnin = 25000,
                            MCMC.samplesize = 5000,
                            parallel = 8,
                            parallel.type = "PSOCK"
                          )
)





# Estrategia incremental
# ===============================================

model_base <- ergm(my_layers ~ 
                     # Solo efectos básicos de densidad y reciprocidad
                     L(~edges + mutual, ~cooperation) +
                     L(~edges + mutual, ~trust) +
                     L(~edges + mutual, ~resources) +
                     L(~edges + mutual, ~values) +
                     L(~edges + mutual, ~kinship),
                   
                   control = control.ergm(
                     MCMLE.maxit = 5,
                     MCMC.samplesize = 1000,
                     MCMLE.termination = "Hummel"
                   )
)


model_homophily <- ergm(my_layers ~ 
                          # Mantenemos efectos base
                          L(~edges + mutual, ~cooperation) +
                          L(~edges + mutual, ~trust) +
                          L(~edges + mutual, ~resources) +
                          L(~edges + mutual, ~values) +
                          L(~edges + mutual, ~kinship) +
                          
                          # Agregamos homofilia
                          L(~nodematch("tipo_org"), ~cooperation) +
                          L(~nodematch("orientacion"), ~cooperation),
                        
                        control = control.ergm(
                          # Usamos los coeficientes del modelo anterior como valores iniciales
                          init = c(coef(model_base), rep(0, 2)),
                          MCMLE.maxit = 15,
                          MCMC.samplesize = 2000,
                          MCMLE.termination = "Hummel"
                        )
)

model_interdep <- ergm(my_layers ~ 
                         # Mantenemos efectos previos
                         L(~edges + mutual, ~cooperation) +
                         L(~edges + mutual, ~trust) +
                         L(~edges + mutual, ~resources) +
                         L(~edges + mutual, ~values) +
                         L(~edges + mutual, ~kinship) +
                         
                         L(~nodematch("tipo_org"), ~cooperation) +
                         L(~nodematch("orientacion"), ~cooperation) +
                         
                         # Agregamos interdependencias principales
                         L(~edges, ~cooperation&trust) +
                         L(~edges, ~cooperation&resources),
                       
                       control = control.ergm(
                         init = c(coef(model_homophily), rep(0, 2)),
                         MCMLE.maxit = 20,
                         MCMC.samplesize = 3000,
                         MCMLE.termination = "Hummel"
                       )
)
summary(model_interdep)


model_final <- ergm(my_layers ~ 
                      # Mantenemos todo lo anterior
                      L(~edges + mutual, ~cooperation) +
                      L(~edges + mutual, ~trust) +
                      L(~edges + mutual, ~resources) +
                      L(~edges + mutual, ~values) +
                      L(~edges + mutual, ~kinship) +
                      
                      L(~nodematch("tipo_org"), ~cooperation) +
                      L(~nodematch("orientacion"), ~cooperation) +
                      
                      # Todas las interdependencias que resultaron significativas
                      L(~edges, ~cooperation&trust) +
                      L(~edges, ~cooperation&resources) +
                      L(~edges, ~trust&resources) +
                      L(~edges, ~cooperation&kinship),
                    
                    control = control.ergm(
                      init = c(coef(model_interdep), rep(0, 2)),  # Ajusta el número de ceros según sea necesario
                      MCMLE.maxit = 30,
                      MCMC.burnin = 75000,
                      MCMC.interval = 1500,
                      MCMC.samplesize = 12500,
                      MCMLE.density.guard = 40,
                      MCMLE.sequential = TRUE,
                      MCMLE.steplength = 1,
                      MCMLE.steplength.margin = 0.5,
                      MCMLE.termination = "Hummel",
                      parallel = 8,
                      parallel.type = "PSOCK"
                    )
)

summary(model_final)

################################################################################
# Latex table for model_final
################################################################################

models <- list(model_final)

# Generate the table in LaTeX format
texreg(models, 
       custom.coef.names = c("L(cooperation) edges", "L(cooperation) mutual", 
                             "L(trust) edges", "L(trust) mutual", 
                             "L(resources) edges", "L(resources) mutual", 
                             "L(values) edges", "L(values) mutual", 
                             "L(kinship) edges", "L(kinship) mutual", 
                             "L(cooperation) homophily.org_type", 
                             "L(cooperation) homophily.orientation", 
                             "L(cooperation - trust) edges", 
                             "L(cooperation - resources) edges", 
                             "L(trust - resources) edges", 
                             "L(cooperation - kinship) edges"),
       file = "output/ergm_results.tex",  # Output file name
       label = "tab:ergm_results",  # Label for referencing in the text
       single.row = TRUE,  # Single row table
       dcolumn = TRUE,  # Column alignment
       use.packages = FALSE,  # Avoid loading packages automatically
       caption = "ERGM results for homophily and interdependences",
       custom.model.names = c("Model 1"),
       override.se = TRUE,  # Override standard errors
       override.coef = TRUE)  # Override coefficients


# Tabla completa 
texreg(list(model_base, model_homophily, model_interdep, model_final),
       custom.coef.names = c(
         "Cooperation: Density",
         "Cooperation: Reciprocity", 
         "Trust: Density",
         "Trust: Reciprocity",
         "Resources: Density", 
         "Resources: Reciprocity",
         "Values: Density",
         "Values: Reciprocity",
         "Kinship: Density",
         "Kinship: Reciprocity",
         "Org. Type Homophily",
         "Orientation Homophily",
         "Coop-Trust Overlap",
         "Coop-Resources Overlap", 
         "Trust-Resources Overlap",
         "Coop-Kinship Overlap"
       ),
       custom.model.names = c("Base", "Homophily", "Core", "Full"),
       stars = c(0.001, 0.01, 0.05),
       digits = 2,
       leading.zero = TRUE,
       caption = "Comparative Analysis of ERGM Models: From Basic to Complex Specifications",
       label = "tab:ergm_comparative",
       booktabs = TRUE,
       dcolumn = TRUE,
       use.packages = TRUE,
       include.rs = FALSE,
       custom.note = "***p < 0.001; **p < 0.01; *p < 0.05",
       custom.additional = list(
         "AIC" = c(4685.20, 4657.48, 3179.39, 2990.31),
         "BIC" = c(4766.12, 4754.58, 3292.68, 3119.78),
         "Log Likelihood" = c(-2332.60, -2316.74, -1575.70, -1479.15)
       ),
       table.placement = "!htb",
       fontsize = "footnotesize",  # Reduce tamaño general de la fuente
       single.row = TRUE          # Coeficientes y SE en una sola fila
       #ci.force = TRUE,            # Fuerza el formato de una sola línea
       #ci.test = 0,                # Elimina los valores p separados
       #ci.force.level = 0.95       # Nivel de confianza para el formato
)



model_final_enhanced <- ergm(my_layers ~ 
                               # Efectos base y existentes (mantenemos igual)
                               L(~edges + mutual, ~cooperation) +
                               L(~edges + mutual, ~trust) +
                               L(~edges + mutual, ~resources) +
                               L(~edges + mutual, ~values) +
                               L(~edges + mutual, ~kinship) +
                               
                               L(~nodematch("tipo_org"), ~cooperation) +
                               L(~nodematch("orientacion"), ~cooperation) +
                               
                               L(~edges, ~cooperation&trust) +
                               L(~edges, ~cooperation&resources) +
                               L(~edges, ~trust&resources) +
                               L(~edges, ~cooperation&kinship) +
                               
                               # Nuevos efectos
                               L(~dgwesp(0.5, fixed=TRUE), ~cooperation) + # Fijamos el parámetro decay
                               L(~dgwesp(0.5, fixed=TRUE), ~resources) +   # Fijamos el parámetro decay
                               L(~nodematch("orientacion"), ~values) +
                               L(~edges, ~values&trust),
                             
                             control = control.ergm(
                               init = c(coef(model_interdep), rep(0, 6)), # Ahora solo 5 nuevos parámetros al fijar decay
                               MCMLE.maxit = 30,
                               MCMC.burnin = 100000,
                               MCMC.interval = 2000,
                               MCMC.samplesize = 15000,
                               MCMLE.density.guard = 40,
                               MCMLE.sequential = TRUE,
                               MCMLE.steplength = 0.5,  # Reducido para mayor estabilidad
                               MCMLE.steplength.margin = 0.25,  # Reducido para mayor estabilidad
                               #MCMLE.termination = "Hummel",
                               parallel = 8,
                               parallel.type = "PSOCK"
                             )
)


