# El presente codigo tiene como objetivo el programar el algoritmo de fuzzy c-means

#' @title Fuzzy C-Means fitness
#' @description esta funcion tiene como objetivo el calcular el fitness de la solucion dada
#' @param U es la matriz de pertenencia
#' @param V es la matriz de centroides
#' @param X es la matriz de datos
#' @param m es el parametro de difuminacion
#' @return el fitness de la solucion dada


fitness_fcm <- function(U, V, X, m) {
  # Calcular la distancia entre los puntos y los centroides
  dist <- as.matrix(dist(rbind(X, V), method = "euclidean"))[1:nrow(X), (nrow(X) + 1):(nrow(X) + nrow(V))]
  
  # Calcular el fitness
  fitness <- sum((U^m) * (dist^2))
  
  return(fitness)
}


#' @title Fuzzy C-Means clustering descenso del gradiente 
#' @description esta funcion tiene como objetivo encontrar una solucion a fuzzy c-means clustering usando descenso del gradiente
#' @param X es la matriz de datos
#' @param k es el numero de clusters
#' @param m es el parametro de difuminacion
#' @param max_iter es el numero maximo de iteraciones
#' @param tol es la tolerancia para la convergencia
#' @return una lista con la matriz de pertenencia, la matriz de centroides y el fitness final


fcm_gradient_descent <- function(X, k, m = 2, max_iter = 100, tol = 1e-5) {
  n <- nrow(X)
  d <- ncol(X)
  
  # Inicializar la matriz de pertenencia U y los centroides V
  U <- matrix(runif(n * k), nrow = n, ncol = k)
  U <- U / rowSums(U)  # Normalizar filas
  
  V <- matrix(0, nrow = k, ncol = d)
  
  for (iter in 1:max_iter) {
    # Actualizar los centroides
    for (j in 1:k) {
      V[j, ] <- colSums((U[, j]^m) * X) / sum(U[, j]^m)
    }
    
    # Actualizar la matriz de pertenencia
    dist <- as.matrix(dist(rbind(X, V), method = "euclidean"))[1:n, (n + 1):(n + k)]
    U_new <- 1 /((dist^2 + .Machine$double.eps)^(2/(m - 1)))
    U_new <- U_new / rowSums(U_new)
    
    # Calcular el cambio en la matriz de pertenencia
    change <- max(abs(U - U_new))
    
    # Actualizar U
    U <- U_new
    
    # Verificar convergencia
    if (change < tol) {
      break
    }
  }
  
  fitness <- fitness_fcm(U, V, X, m)
  
  return(list(U = U, V = V, fitness = fitness))
}


#' @titlle cruzamiento de la poblacion
#' @description esta funcion tiene como objetivo cruzar la poblacion de soluciones
#' @param X es la matriz de datos
#' @param population es la poblacion de soluciones
#' @param probabilities es el vector de probabilidades de cruce
#' @param k es el numero de clusters
#' @param m es el parametro de difuminacion
#' @return una lista con la poblacion cruzada
cruce_population <- function(X, population, probabilities, k, m, best_population) {
  pop_size <- length(population)
  new_population <- vector("list", pop_size + 2)
  
  for (j in seq(1, pop_size, by = 2)) {
    
    # Seleccionar dos padres basados en las probabilidades
    parents_indices <- sample(1:pop_size, size = 2, prob = probabilities)
    parent1 <- population[[parents_indices[1]]]
    parent2 <- population[[parents_indices[2]]]
    
    # Realizar el cruce
    V1 <- parent1$V*1000
    V2 <- parent2$V*1000
    
    #Se crean los hijos
    V_h1 <- matrix(0, nrow = nrow(V1), ncol = ncol(V1))
    V_h2 <- matrix(0, nrow = nrow(V1), ncol = ncol(V1))

    for(i in 1:k){
        # Convertir los centroides a bits
        centroP1 <- intToBits(V1[i, ])

        
        centroP2 <- intToBits(V2[i, ])

        #Se procede a realizar el cruce
        centroHijo1 <- intToBits(0)
        centroHijo2 <- intToBits(0)

        cruce <- sample(c(TRUE, FALSE), size = length(centroP1), replace = TRUE)
        centroHijo1[cruce] <- centroP1[cruce]
        centroHijo1[!cruce] <- centroP2[!cruce]

        V_h1[i,] <- packBits(centroHijo1, type = "integer")/1000


        
        centroHijo2[cruce] <- centroP2[cruce]
        centroHijo2[!cruce] <- centroP1[!cruce]

        V_h2[i,] <- packBits(centroHijo2, type = "integer")/1000
    }
    if (j + 1 > pop_size) {
        #Se incluye el primer hijo en la nueva poblacion
        #Se calculan los nuevos grados de pertenencia para el primer hijo
        dist <- as.matrix(dist(rbind(X, V_h1), method = "euclidean"))[1:nrow(X), (nrow(X) + 1):(nrow(X) + nrow(V_h1))]
        U_new <- 1 / ((dist^2 + .Machine$double.eps)^(2/(m - 1)))
        U_new <- U_new / rowSums(U_new)
        
        # Crear la nueva solucion
        new_population[[j]] <- list(U = U_new, V = V_h1, fitness = fitness_fcm(U_new, V_h1, X, m))

    }else{
        #Se calculan los nuevos grados de pertenencia para el primer hijo
        dist <- as.matrix(dist(rbind(X, V_h1), method = "euclidean"))[1:nrow(X), (nrow(X) + 1):(nrow(X) + nrow(V_h1))]
        U_new <- 1 / ((dist^2 + .Machine$double.eps)^(2/(m - 1)))
        U_new <- U_new / rowSums(U_new)
        
        # Crear la nueva solucion
        new_population[[j]] <- list(U = U_new, V = V_h1, fitness = fitness_fcm(U_new, V_h1, X, m))
        
        #Se calculan los nuevos grados de pertenencia para el segundo hijo
        dist <- as.matrix(dist(rbind(X, V_h2), method = "euclidean"))[1:nrow(X), (nrow(X) + 1):(nrow(X) + nrow(V_h2))]
        U_new <- 1 / ((dist^2 + .Machine$double.eps)^(2/(m - 1)))
        U_new <- U_new / rowSums(U_new)
        
        # Crear la nueva solucion
        new_population[[j+1]] <- list(U = U_new, V = V_h2, fitness = fitness_fcm(U_new, V_h2, X, m))
    }
  }
  new_population[[length(new_population)]] <- best_population[[1]]
  new_population[[length(new_population) - 1]] <- best_population[[2]]
  return(new_population)
}




#' @title fuzzy c-means clustering algoritmo genetico
#' @description esta funcion tiene como objetivo encontrar una solucion a fuzzy c-means clustering usando
#' @param X es la matriz de datos
#' @param k es el numero de clusters
#' @param m es el parametro de difuminacion
#' @param pop_size es el tamaño de la poblacion
#' @param max_iter es el numero maximo de iteraciones
#' @param tol es la tolerancia para la convergencia
#' @return una lista con la mejor solucion encontrada, la matriz de pertenencia, la matriz de centroides y el fitness final

fcm_ga <- function(X, k, m = 2, pop_size = 50, max_iter = 100, tol = 1e-5, setps = 10) {
  n <- nrow(X)

  # Inicializar la poblacion
  population <- vector("list", pop_size)
  for (i in 1:pop_size) {
    population[[i]] <- fcm_gradient_descent(X, k, m, max_iter, tol)
  }
  for(i in 1:max_iter){
    
    # Evaluar la poblacion
    fitness_values <- sapply(population, function(sol) sol$fitness)
    
    # Seleccionar los mejores individuos
    best_indices <- order(fitness_values)[1:2]
    best_population <- population[best_indices]
    if(i > setps){
        if(abs(best_solution$fitness - best_population[[1]]$fitness)< tol) {
        # Si la mejor solucion es suficientemente buena, detener el algoritmo
        break
        }
    }
    
    
    # Actualizar la mejor solucion
    best_solution <- best_population[[1]]
    
    # Calcular las probabilidades de cruce
    probabilities <- fitness_values / sum(fitness_values)
    
    # Cruzar la poblacion
    population <- cruce_population(X, population, probabilities, k, m, best_population)
  }
  
  return(best_solution)
}
