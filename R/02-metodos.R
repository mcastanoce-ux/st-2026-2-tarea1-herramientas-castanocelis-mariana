
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

# Validacion comun de la serie de entrada.
.validar_y <- function(y, minimo = 3L) {
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y) && is.null(dim(y)),
    "y no puede contener NA: la serie debe venir sin huecos" = !any(is.na(y)),
    "y debe contener valores finitos" = all(is.finite(y))
  )
  if (length(y) < minimo) {
    stop(sprintf("la serie tiene %d observaciones y se requieren al menos %d.",
                 length(y), minimo), call. = FALSE)
  }
  invisible(TRUE)
}

.validar_constante <- function(a, nombre) {
  if (!(is.numeric(a) && length(a) == 1L && is.finite(a) && a > 0 && a < 1)) {
    stop(sprintf("la constante %s debe ser un numero en el intervalo abierto (0, 1); se recibio %s.",
                 nombre, paste(a, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

.validar_ventana <- function(k, n) {
  if (!(is.numeric(k) && length(k) == 1L && k == as.integer(k) && k >= 2)) {
    stop("la ventana k debe ser un entero mayor o igual que 2.", call. = FALSE)
  }
  if (k > n) {
    stop(sprintf("la ventana k = %d es mayor que la longitud de la serie (%d).",
                 as.integer(k), n), call. = FALSE)
  }
  invisible(TRUE)
}

# --------------------------------------------------------------------------
# (a) Metodos de nivel
# --------------------------------------------------------------------------

ajustar_media <- function(y) {
  .validar_y(y, minimo = 3L)
  T_n <- length(y)
  yhat <- rep(NA_real_, T_n)
  suma <- 0
  for (t in seq_len(T_n - 1L)) {
    suma <- suma + y[t]               # actualizacion recursiva de la suma
    yhat[t + 1L] <- suma / t          # media de Y_1..Y_t
  }
  suma_total <- suma + y[T_n]
  media_final <- suma_total / T_n
  list(
    yhat = yhat,
    pronosticar = function(h) {
      stopifnot("h debe ser un entero positivo" = is.numeric(h) && h >= 1)
      rep(media_final, as.integer(h))
    },
    parametros = list(
      metodo = "media simple",
      calentamiento = 1L,
      p = 1L,                          # un parametro estimado: el nivel
      media_final = media_final,
      T = T_n
    )
  )
}

#Media movil
ajustar_mm <- function(y, k) {
  .validar_y(y, minimo = 3L)
  .validar_ventana(k, length(y))
  k <- as.integer(k)
  T_n <- length(y)
  yhat <- rep(NA_real_, T_n)
  mm   <- rep(NA_real_, T_n)          # MM_t(k)
  suma <- 0
  for (t in seq_len(T_n)) {
    suma <- suma + y[t]
    if (t > k) suma <- suma - y[t - k]     # ventana movil: un solo recorrido
    if (t >= k) {
      mm[t] <- suma / k
      if (t < T_n) yhat[t + 1L] <- mm[t]
    }
  }
  mm_final <- mm[T_n]
  list(
    yhat = yhat,
    pronosticar = function(h) {
      stopifnot("h debe ser un entero positivo" = is.numeric(h) && h >= 1)
      rep(mm_final, as.integer(h))
    },
    parametros = list(
      metodo = "media movil",
      k = k,
      calentamiento = k,
      p = 1L,
      mm = mm,
      mm_final = mm_final,
      T = T_n
    )
  )
}

#' Suavizamiento exponencial simple, escrito en forma de correccion de error
ajustar_ses <- function(y, alpha) {
  .validar_y(y, minimo = 3L)
  .validar_constante(alpha, "alpha")
  T_n <- length(y)
  yhat <- rep(NA_real_, T_n)
  yhat[2L] <- y[1L]
  for (t in 2:(T_n - 1L)) {
    e <- y[t] - yhat[t]
    yhat[t + 1L] <- yhat[t] + alpha * e            # forma de correccion de error
  }
  e_T <- y[T_n] - yhat[T_n]
  nivel_final <- yhat[T_n] + alpha * e_T           # = alpha Y_T + (1-alpha) Yhat_T
  
  # Verificacion interna: forma de correccion de error vs promedio ponderado
  yhat_pp <- rep(NA_real_, T_n)
  yhat_pp[2L] <- y[1L]
  for (t in 2:(T_n - 1L)) {
    yhat_pp[t + 1L] <- alpha * y[t] + (1 - alpha) * yhat_pp[t]
  }
  dif_max <- max(abs(yhat - yhat_pp), na.rm = TRUE)
  if (dif_max > 1e-10) {
    stop("ajustar_ses(): la forma de correccion de error no coincide con el promedio ponderado.",
         call. = FALSE)
  }
  
  list(
    yhat = yhat,
    pronosticar = function(h) {
      stopifnot("h debe ser un entero positivo" = is.numeric(h) && h >= 1)
      rep(nivel_final, as.integer(h))
    },
    parametros = list(
      metodo = "suavizamiento exponencial simple",
      alpha = alpha,
      calentamiento = 1L,
      p = 1L,
      nivel_final = nivel_final,
      dif_max_formas = dif_max,
      T = T_n
    )
  )
}

#' Doble media movil de orden k.
ajustar_dmm <- function(y, k) {
  .validar_y(y, minimo = 3L)
  .validar_ventana(k, length(y))
  k <- as.integer(k)
  T_n <- length(y)
  if (2L * k - 1L >= T_n) {
    stop(sprintf("la ventana k = %d exige %d periodos de calentamiento y la serie tiene %d observaciones.",
                 k, 2L * k - 1L, T_n), call. = FALSE)
  }
  mm  <- rep(NA_real_, T_n)
  dmm <- rep(NA_real_, T_n)
  E   <- rep(NA_real_, T_n)
  b1  <- rep(NA_real_, T_n)
  yhat <- rep(NA_real_, T_n)
  
  suma1 <- 0; suma2 <- 0
  for (t in seq_len(T_n)) {
    suma1 <- suma1 + y[t]
    if (t > k) suma1 <- suma1 - y[t - k]
    if (t >= k) {
      mm[t] <- suma1 / k
      suma2 <- suma2 + mm[t]
      if (t >= 2L * k) suma2 <- suma2 - mm[t - k]
      if (t >= 2L * k - 1L) {
        dmm[t] <- suma2 / k
        E[t]  <- 2 * mm[t] - dmm[t]
        b1[t] <- (2 / (k - 1)) * (mm[t] - dmm[t])
        if (t < T_n) yhat[t + 1L] <- E[t] + b1[t]
      }
    }
  }
  E_T <- E[T_n]; b1_T <- b1[T_n]
  list(
    yhat = yhat,
    pronosticar = function(h) {
      stopifnot("h debe ser un entero positivo" = is.numeric(h) && h >= 1)
      E_T + b1_T * seq_len(as.integer(h))
    },
    parametros = list(
      metodo = "doble media movil",
      k = k,
      calentamiento = 2L * k - 1L,
      p = 2L,                       # nivel y pendiente
      mm = mm, dmm = dmm, E = E, b1 = b1,
      E_final = E_T, b1_final = b1_T,
      T = T_n
    )
  )
}

# --------------------------------------------------------------------------
# (b) Tendencias por minimos cuadrados
# --------------------------------------------------------------------------

# Errores estandar HAC con nucleo de Bartlett y L = floor(4 (T/100)^(2/9)).
.hac_bartlett <- function(X, u) {
  T_n <- nrow(X); kx <- ncol(X)
  L <- floor(4 * (T_n / 100)^(2 / 9))
  XtX_inv <- solve(crossprod(X))
  Omega <- matrix(0, kx, kx)
  for (t in seq_len(T_n)) {
    xt <- matrix(X[t, ], ncol = 1)
    Omega <- Omega + (u[t]^2) * (xt %*% t(xt))
  }
  if (L >= 1) {
    for (l in seq_len(L)) {
      w <- 1 - l / (L + 1)                       # peso de Bartlett
      Gl <- matrix(0, kx, kx)
      for (t in (l + 1):T_n) {
        xt  <- matrix(X[t, ], ncol = 1)
        xtl <- matrix(X[t - l, ], ncol = 1)
        Gl <- Gl + u[t] * u[t - l] * (xt %*% t(xtl) + xtl %*% t(xt))
      }
      Omega <- Omega + w * Gl
    }
  }
  V <- XtX_inv %*% Omega %*% XtX_inv
  list(V = V, ee = sqrt(pmax(diag(V), 0)), rezagos = L)
}


ajustar_tendencia <- function(y, tipo = c("lineal", "cuadratica", "exponencial"),
                              corregir_sesgo = FALSE) {
  tipo <- match.arg(tipo)
  .validar_y(y, minimo = 5L)
  stopifnot("corregir_sesgo debe ser TRUE o FALSE" = is.logical(corregir_sesgo) &&
              length(corregir_sesgo) == 1L)
  T_n <- length(y)
  t_idx <- seq_len(T_n)
  
  if (tipo == "exponencial" && any(y <= 0)) {
    stop("ajustar_tendencia(): la tendencia exponencial se estima sobre log(Y_t) y la serie tiene valores no positivos.",
         call. = FALSE)
  }
  
  # Matriz de diseno construida a mano
  X <- switch(tipo,
              lineal      = cbind(1, t_idx),
              cuadratica  = cbind(1, t_idx, t_idx^2),
              exponencial = cbind(1, t_idx)
  )
  colnames(X) <- switch(tipo,
                        lineal      = c("b0", "b1"),
                        cuadratica  = c("b0", "b1", "b2"),
                        exponencial = c("a", "theta")
  )
  z <- if (tipo == "exponencial") log(y) else y
  
  beta <- solve(crossprod(X), crossprod(X, z))     # ecuaciones normales
  beta <- as.numeric(beta)
  names(beta) <- colnames(X)
  
  z_ajustado <- as.numeric(X %*% beta)
  u <- z - z_ajustado                              # residuos en la escala estimada
  kx <- ncol(X)
  gl <- T_n - kx
  sigma2 <- sum(u^2) / gl
  R2 <- 1 - sum(u^2) / sum((z - mean(z))^2)
  
  ee_ols <- sqrt(diag(solve(crossprod(X))) * sigma2)
  hac <- .hac_bartlett(X, u)
  t_ols <- beta / ee_ols
  t_hac <- beta / hac$ee
  tabla <- tibble::tibble(
    coeficiente = names(beta),
    estimacion  = beta,
    ee_ordinario = ee_ols,
    t_ordinario  = t_ols,
    p_ordinario  = 2 * stats::pt(abs(t_ols), df = gl, lower.tail = FALSE),
    ee_hac       = hac$ee,
    t_hac        = t_hac,
    p_hac        = 2 * stats::pt(abs(t_hac), df = gl, lower.tail = FALSE)
  )
  dw <- sum(diff(u)^2) / sum(u^2)
  
  factor_sesgo <- if (tipo == "exponencial" && corregir_sesgo) exp(sigma2 / 2) else 1
  yhat <- if (tipo == "exponencial") exp(z_ajustado) * factor_sesgo else z_ajustado
  
  coef <- beta
  list(
    yhat = yhat,
    pronosticar = function(h) {
      stopifnot("h debe ser un entero positivo" = is.numeric(h) && h >= 1)
      hh <- as.integer(h)
      tf <- T_n + seq_len(hh)
      switch(tipo,
             lineal      = coef[1] + coef[2] * tf,
             cuadratica  = coef[1] + coef[2] * tf + coef[3] * tf^2,
             exponencial = exp(coef[1] + coef[2] * tf) * factor_sesgo
      )
    },
    parametros = list(
      metodo = paste("tendencia", tipo),
      tipo = tipo,
      calentamiento = 0L,
      p = kx,
      coeficientes = coef,
      b0_exponencial = if (tipo == "exponencial") exp(coef[1]) else NA_real_,
      b1_exponencial = if (tipo == "exponencial") exp(coef[2]) else NA_real_,
      corregir_sesgo = corregir_sesgo,
      factor_sesgo = factor_sesgo,
      tabla_coeficientes = tabla,
      R2 = R2,
      sigma2 = sigma2,
      gl = gl,
      rezagos_hac = hac$rezagos,
      durbin_watson = dw,
      residuos_escala_estimacion = u,
      T = T_n
    )
  )
}

# --------------------------------------------------------------------------
# (c) Holt lineal
# --------------------------------------------------------------------------

#' Holt lineal con L_1 = Y_1 y That_1 = 0
ajustar_holt <- function(y, alpha, beta) {
  .validar_y(y, minimo = 4L)
  .validar_constante(alpha, "alpha")
  .validar_constante(beta, "beta")
  T_n <- length(y)
  L <- rep(NA_real_, T_n); Tt <- rep(NA_real_, T_n); yhat <- rep(NA_real_, T_n)
  L[1L] <- y[1L]; Tt[1L] <- 0
  for (t in 2:T_n) {
    yhat[t] <- L[t - 1L] + Tt[t - 1L]
    e <- y[t] - yhat[t]
    L[t]  <- alpha * y[t] + (1 - alpha) * yhat[t]
    Tt[t] <- beta * (L[t] - L[t - 1L]) + (1 - beta) * Tt[t - 1L]
  }
  # Verificacion de la forma de correccion de error
  L2 <- rep(NA_real_, T_n); T2 <- rep(NA_real_, T_n)
  L2[1L] <- y[1L]; T2[1L] <- 0
  for (t in 2:T_n) {
    e <- y[t] - (L2[t - 1L] + T2[t - 1L])
    L2[t] <- L2[t - 1L] + T2[t - 1L] + alpha * e
    T2[t] <- T2[t - 1L] + alpha * beta * e
  }
  dif_max <- max(c(abs(L - L2), abs(Tt - T2)), na.rm = TRUE)
  if (dif_max > 1e-10) {
    stop("ajustar_holt(): la forma de correccion de error no coincide con las ecuaciones de nivel y pendiente.",
         call. = FALSE)
  }
  L_T <- L[T_n]; T_T <- Tt[T_n]
  list(
    yhat = yhat,
    pronosticar = function(h) {
      stopifnot("h debe ser un entero positivo" = is.numeric(h) && h >= 1)
      L_T + T_T * seq_len(as.integer(h))
    },
    parametros = list(
      metodo = "Holt lineal",
      alpha = alpha, beta = beta,
      calentamiento = 1L,
      p = 2L,
      nivel = L, pendiente = Tt,
      nivel_final = L_T, pendiente_final = T_T,
      dif_max_correccion_error = dif_max,
      T = T_n
    )
  )
}

# --------------------------------------------------------------------------
# (d) Optimizacion de constantes y ventanas
# --------------------------------------------------------------------------

optimizar <- function(y, metodo = c("mm", "dmm", "ses", "holt"), rejilla = NULL) {
  metodo <- match.arg(metodo)
  .validar_y(y, minimo = 5L)
  if (is.null(rejilla)) {
    rejilla <- switch(metodo,
                      mm   = 2:12,
                      dmm  = 2:12,
                      ses  = seq(0.02, 0.98, by = 0.02),
                      holt = expand.grid(alpha = seq(0.05, 0.95, by = 0.05),
                                         beta  = seq(0.05, 0.95, by = 0.05))
    )
  }
  mse_de <- function(ajuste) {
    e <- y - ajuste$yhat
    mean(e^2, na.rm = TRUE)
  }
  if (metodo == "holt") {
    stopifnot("la rejilla de Holt debe tener columnas alpha y beta" =
                all(c("alpha", "beta") %in% names(rejilla)))
    mse <- vapply(seq_len(nrow(rejilla)), function(i) {
      tryCatch(mse_de(ajustar_holt(y, rejilla$alpha[i], rejilla$beta[i])),
               error = function(e) NA_real_)
    }, numeric(1))
    tabla <- tibble::as_tibble(rejilla) %>% dplyr::mutate(mse = mse)
    i <- which.min(tabla$mse)
    optimo <- list(alpha = tabla$alpha[i], beta = tabla$beta[i], mse = tabla$mse[i])
    borde <- tabla$alpha[i] %in% range(rejilla$alpha) || tabla$beta[i] %in% range(rejilla$beta)
  } else {
    valores <- as.numeric(rejilla)
    mse <- vapply(valores, function(v) {
      tryCatch(
        mse_de(switch(metodo,
                      mm  = ajustar_mm(y, v),
                      dmm = ajustar_dmm(y, v),
                      ses = ajustar_ses(y, v))),
        error = function(e) NA_real_)
    }, numeric(1))
    nombre <- if (metodo == "ses") "alpha" else "k"
    tabla <- tibble::tibble(valor = valores, mse = mse)
    names(tabla)[1] <- nombre
    i <- which.min(tabla$mse)
    optimo <- list(valor = valores[i], mse = tabla$mse[i])
    names(optimo)[1] <- nombre
    borde <- valores[i] %in% range(valores)
  }
  list(
    metodo = metodo,
    rejilla = tabla,
    optimo = optimo,
    en_borde = borde
  )
}

#' Pronostico ingenuo (referente obligatorio)
ajustar_ingenuo <- function(y, s = 1) {
  .validar_y(y, minimo = 2L)
  stopifnot("s debe ser un entero positivo menor que la longitud de y" =
              is.numeric(s) && s >= 1 && s < length(y))
  s <- as.integer(s)
  T_n <- length(y)
  yhat <- rep(NA_real_, T_n)
  if (s < T_n) yhat[(s + 1L):T_n] <- y[1:(T_n - s)]
  list(
    yhat = yhat,
    pronosticar = function(h) {
      hh <- as.integer(h)
      if (s == 1L) return(rep(y[T_n], hh))
      idx <- T_n + seq_len(hh) - s * (floor((seq_len(hh) - 1) / s) + 1)
      y[idx]
    },
    parametros = list(metodo = if (s == 1L) "ingenuo" else "ingenuo estacional",
                      s = s, p = 0L, calentamiento = s, T = T_n)
  )
}
