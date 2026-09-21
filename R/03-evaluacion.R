suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})


medidas <- function(y, yhat, escala = NULL) {
  stopifnot(
    "y y yhat deben tener la misma longitud" = length(y) == length(yhat),
    "y debe ser numerico" = is.numeric(y), "yhat debe ser numerico" = is.numeric(yhat)
  )
  ok <- !is.na(yhat) & !is.na(y)
  e <- y[ok] - yhat[ok]
  n <- length(e)
  stopifnot("no hay pares (y, yhat) utilizables" = n >= 1)
  mape <- if (any(y[ok] == 0)) NA_real_ else mean(abs(e / y[ok])) * 100
  tibble::tibble(
    n    = n,
    MSE  = mean(e^2),
    RMSE = sqrt(mean(e^2)),
    MAD  = mean(abs(e)),
    MAPE = mape,
    MASE = if (is.null(escala) || !is.finite(escala) || escala <= 0) NA_real_
    else mean(abs(e)) / escala
  )
}

#' Escalador del MASE: MAE de un paso del ingenuo dentro del tramo de estimacion.

escala_mase <- function(y_estimacion, s = 1) {
  stopifnot("y_estimacion debe ser numerico" = is.numeric(y_estimacion))
  s <- as.integer(s)
  stopifnot("s debe ser menor que la longitud del tramo" = s < length(y_estimacion))
  mean(abs(y_estimacion[(s + 1):length(y_estimacion)] -
             y_estimacion[1:(length(y_estimacion) - s)]))
}

#Prueba de Ljung-Box.
ljung_box <- function(r, T, m = length(r), p = 0) {
  stopifnot(
    "r debe ser numerico" = is.numeric(r),
    "T debe ser un entero positivo" = is.numeric(T) && T > 1,
    "m no puede exceder la longitud de r" = m <= length(r),
    "p debe ser un entero no negativo" = p >= 0,
    "m debe ser mayor que p para tener grados de libertad positivos" = (m - p) > 0
  )
  h <- seq_len(m)
  Q <- T * (T + 2) * sum(r[h]^2 / (T - h))
  gl <- m - p
  list(
    prueba = "Ljung-Box",
    estadistico = Q,
    gl = gl,
    critico = stats::qchisq(0.95, df = gl),
    valor_p = stats::pchisq(Q, df = gl, lower.tail = FALSE),
    H0 = "rho_1 = ... = rho_m = 0 (los errores/observaciones no estan autocorrelacionados)",
    H1 = "algun rho_h distinto de cero para h = 1,...,m",
    formula = "Q_m = T(T+2) sum_{h=1}^{m} r_h^2/(T-h) ~ chi^2_{m-p}",
    region = sprintf("se rechaza H0 si Q_m > chi^2_{0.95, %d} = %.4f", gl,
                     stats::qchisq(0.95, df = gl))
  )
}

#' Prueba de Jarque-Bera. JB = (N/6)(A^2 + (K-3)^2/4) ~ chi^2_2.
jarque_bera <- function(e) {
  e <- e[!is.na(e)]
  stopifnot("e debe ser numerico" = is.numeric(e),
            "se requieren al menos 4 errores" = length(e) >= 4L)
  N <- length(e)
  d <- e - mean(e)
  m2 <- sum(d^2) / N; m3 <- sum(d^3) / N; m4 <- sum(d^4) / N
  A <- m3 / m2^(3 / 2)          # asimetria
  K <- m4 / m2^2                # curtosis (3 bajo normalidad)
  JB <- (N / 6) * (A^2 + (K - 3)^2 / 4)
  list(
    prueba = "Jarque-Bera",
    estadistico = JB,
    gl = 2,
    critico = stats::qchisq(0.95, df = 2),
    valor_p = stats::pchisq(JB, df = 2, lower.tail = FALSE),
    asimetria = A, curtosis = K, N = N,
    H0 = "los errores provienen de una distribucion normal (A = 0 y K = 3)",
    H1 = "los errores no son normales (A != 0 o K != 3)",
    formula = "JB = (N/6)(A^2 + (K-3)^2/4) ~ chi^2_2",
    region = sprintf("se rechaza H0 si JB > chi^2_{0.95,2} = %.4f", stats::qchisq(0.95, 2)),
    advertencia = if (N < 20)
      sprintf("JB es asintotica y aqui solo hay %d errores: la decision se lee con reserva.", N)
    else NA_character_
  )
}

#' Estadistico de Durbin-Watson. 
durbin_watson <- function(e) {
  e <- e[!is.na(e)]
  stopifnot("e debe ser numerico" = is.numeric(e),
            "se requieren al menos 3 errores" = length(e) >= 3L)
  d <- sum(diff(e)^2) / sum(e^2)
  list(
    prueba = "Durbin-Watson",
    estadistico = d,
    gl = NA_real_,
    critico = NA_real_,      
    valor_p = NA_real_,
    rho_aprox = 1 - d / 2,
    n = length(e),
    H0 = "rho = 0 (errores no autocorrelacionados de orden 1)",
    H1 = "rho > 0 (autocorrelacion positiva de primer orden)",
    formula = "d = sum_t (e_t - e_{t-1})^2 / sum_t e_t^2, con d aprox 2(1 - rho)",
    region = "se rechaza H0 si d < dL; no se concluye si dL <= d <= dU; no se rechaza si d > dU (cotas de la tabla para T y k regresores)"
  )
}

#Prueba t de media cero sobre los errores.
t_media_cero <- function(e) {
  e <- e[!is.na(e)]
  stopifnot("se requieren al menos 3 errores" = length(e) >= 3L)
  n <- length(e); m <- mean(e); s <- stats::sd(e)
  tt <- m / (s / sqrt(n))
  list(
    prueba = "t de media cero",
    estadistico = tt,
    gl = n - 1,
    critico = stats::qt(0.975, df = n - 1),
    valor_p = 2 * stats::pt(abs(tt), df = n - 1, lower.tail = FALSE),
    media = m, ee = s / sqrt(n), n = n,
    H0 = "mu_e = 0 (el metodo no esta sesgado)",
    H1 = "mu_e != 0 (el metodo esta sesgado)",
    formula = "t = ebar / (s_e/sqrt(n)) ~ t_{n-1}",
    region = sprintf("se rechaza H0 si |t| > t_{0.975, %d} = %.4f", n - 1,
                     stats::qt(0.975, df = n - 1))
  )
}

#Prueba t individual sobre un coeficiente de tendencia con error estandar

t_coeficiente <- function(estimacion, ee, gl, nombre = "beta_j") {
  tt <- estimacion / ee
  list(
    prueba = sprintf("t sobre %s (EE robusto HAC)", nombre),
    estadistico = tt,
    gl = gl,
    critico = stats::qt(0.975, df = gl),
    valor_p = 2 * stats::pt(abs(tt), df = gl, lower.tail = FALSE),
    H0 = sprintf("%s = 0", nombre),
    H1 = sprintf("%s != 0", nombre),
    formula = sprintf("t = %s_hat / EE_HAC(%s_hat) ~ t_{T-k}", nombre, nombre),
    region = sprintf("se rechaza H0 si |t| > t_{0.975, %d} = %.4f", gl,
                     stats::qt(0.975, df = gl))
  )
}

#Interpretacion de las hipotesis
reportar_prueba <- function(prueba, lectura = NULL, alfa = 0.05) {
  stopifnot("prueba debe ser una lista devuelta por una funcion de prueba" =
              is.list(prueba) && !is.null(prueba$estadistico))
  decision <- if (is.na(prueba$valor_p)) {
    "la decision se toma con las cotas dL y dU de la tabla (ver informe)"
  } else if (prueba$valor_p < alfa) {
    sprintf("se rechaza H0 al %g %%", 100 * alfa)
  } else {
    sprintf("no se rechaza H0 al %g %%", 100 * alfa)
  }
  cat("\n--- ", prueba$prueba, " ---\n", sep = "")
  cat("1. H0: ", prueba$H0, "\n   H1: ", prueba$H1, "\n", sep = "")
  cat("2. Estadistico: ", prueba$formula,
      if (!is.na(prueba$gl[1])) sprintf("  (gl = %s)", paste(prueba$gl, collapse = ", ")) else "",
      "\n", sep = "")
  cat("3. Region de rechazo al 5 %: ", prueba$region, "\n", sep = "")
  cat(sprintf("4. Valor observado = %.6f ; valor p = %s\n",
              prueba$estadistico,
              if (is.na(prueba$valor_p)) "no disponible (ver cotas)" else
                formatC(prueba$valor_p, format = "g", digits = 4)))
  cat("5. Decision: ", decision, "\n", sep = "")
  cat("6. Lectura: ", if (is.null(lectura)) "(se escribe en el informe)" else lectura, "\n", sep = "")
  if (!is.null(prueba$advertencia) && !is.na(prueba$advertencia)) {
    cat("   Advertencia: ", prueba$advertencia, "\n", sep = "")
  }
  
  
  x <- prueba
  invisible(tibble::tibble(
    prueba = x$prueba,
    estadistico = x$estadistico,
    gl = paste(x$gl, collapse = ", "),
    critico = x$critico,
    valor_p = x$valor_p,
    decision = decision
  ))
}

#' Validacion completa de los errores de un paso de un metodo.

validar_errores <- function(e, fechas = NULL, p = 0, m = NULL, titulo = "Errores de un paso",
                            imprimir = TRUE) {
  stopifnot("e debe ser numerico" = is.numeric(e))
  ok <- !is.na(e)
  e_ok <- e[ok]
  n <- length(e_ok)
  stopifnot("se requieren al menos 8 errores para validar" = n >= 8L)
  if (is.null(m)) m <- min(floor(n / 4), 24)
  if (m <= p) m <- p + 1L
  
  fechas_ok <- if (is.null(fechas)) seq_len(n) else fechas[ok]
  d <- tibble::tibble(x = fechas_ok, e = e_ok)
  
  g_e <- ggplot(d, aes(x = x, y = e)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_line(colour = "grey25") + geom_point(size = 0.8) +
    labs(title = titulo, x = "Fecha", y = "Error de un paso",
         caption = sprintf("n = %d errores.", n)) +
    theme_bw(base_size = 10) + theme(plot.caption = element_text(hjust = 0, size = 8))
  
  g_corr <- correlograma(e_ok, m = m, titulo = "Correlograma de los errores")
  r <- attr(g_corr, "tabla")$r
  
  pruebas <- list(
    t_media = t_media_cero(e_ok),
    ljung_box = ljung_box(r, T = n, m = m, p = p),
    jarque_bera = jarque_bera(e_ok),
    durbin_watson = durbin_watson(e_ok)
  )
  if (imprimir) {
    cat("\n=== Validacion de errores:", titulo, "===\n")
    cat(sprintf("Numero de errores: %d. Rezagos usados: m = %d. Parametros del metodo: p = %d (gl = %d).\n",
                n, m, p, m - p))
    invisible(lapply(pruebas, reportar_prueba))
  }
  resumen <- dplyr::bind_rows(lapply(pruebas, function(pr) tibble::tibble(
    prueba = pr$prueba, estadistico = pr$estadistico,
    gl = paste(pr$gl, collapse = ","), critico = pr$critico, valor_p = pr$valor_p)))
  
  list(
    n = n, m = m, p = p,
    grafico_errores = g_e,
    correlograma = g_corr,
    banda = attr(g_corr, "banda"),
    pruebas = pruebas,
    resumen = resumen
  )
}

#' Grafico de la curva de MSE (una constante) o del mapa de MSE (alpha, beta).
graficar_optimizacion <- function(opt, titulo = "MSE sobre la rejilla") {
  stopifnot("opt debe venir de optimizar()" = is.list(opt) && !is.null(opt$rejilla))
  tabla <- opt$rejilla
  if (opt$metodo == "holt") {
    ggplot(tabla, aes(x = alpha, y = beta, fill = mse)) +
      geom_tile() +
      geom_point(data = tibble::tibble(alpha = opt$optimo$alpha, beta = opt$optimo$beta),
                 aes(x = alpha, y = beta), inherit.aes = FALSE,
                 shape = 4, size = 3, stroke = 1.2, colour = "white") +
      scale_fill_viridis_c(option = "C") +
      labs(title = titulo, x = expression(alpha), y = expression(beta), fill = "MSE (log)",
           caption = sprintf("Optimo: alpha = %.2f, beta = %.2f, MSE = %.4f.%s",
                             opt$optimo$alpha, opt$optimo$beta, opt$optimo$mse,
                             if (opt$en_borde) " El optimo cae en el borde de la rejilla." else "")) +
      theme_bw(base_size = 10) + theme(plot.caption = element_text(hjust = 0, size = 8))
  } else {
    nombre <- names(tabla)[1]
    ggplot(tabla, aes(x = .data[[nombre]], y = mse)) +
      geom_line(colour = "grey25") + geom_point(size = 1) +
      geom_point(data = tabla[which.min(tabla$mse), ], colour = "#c0392b", size = 3) +
      labs(title = titulo, x = nombre, y = "MSE de un paso (tramo de estimacion)",
           caption = sprintf("Optimo: %s = %s, MSE = %.4f.%s",
                             nombre, format(opt$optimo[[1]]), opt$optimo$mse,
                             if (opt$en_borde) " El optimo cae en el borde de la rejilla." else "")) +
      theme_bw(base_size = 10) + theme(plot.caption = element_text(hjust = 0, size = 8))
  }
}

# ---------------------------------------------------------------------------
# Protocolo de un ejemplo 
# ---------------------------------------------------------------------------


protocolo_ejemplo <- function(datos, nombre, ajustador, etiqueta, p,
                              estacional = FALSE, h = NULL, imprimir = TRUE) {
  stopifnot("datos debe venir de leer_serie()" = all(c("t", "fecha", "y") %in% names(datos)))
  s <- if (estacional) attr(datos, "frecuencia") else 1
  part <- particionar(datos, h = h, estacional = estacional)
  y_est <- part$estimacion$y
  y_val <- part$validacion$y
  hh <- part$h
  
  # 1. Descripcion y correlograma de la serie completa
  descripcion <- describir_serie(datos, nombre)
  g_serie <- graficar_serie(datos, sprintf("%s: serie observada", nombre))
  corr_serie <- correlograma(datos, titulo = sprintf("%s: correlograma de la serie", nombre))
  rh_serie <- contrastar_rh(datos$y)
  m_serie <- nrow(attr(corr_serie, "tabla"))
  lb_serie <- ljung_box(attr(corr_serie, "tabla")$r, T = nrow(datos), m = m_serie, p = 0)
  
  # 2. Ajuste sobre el tramo de estimacion
  ajuste <- ajustador(y_est)
  e_est <- y_est - ajuste$yhat
  pronos <- ajuste$pronosticar(hh)
  
  # 3. Referente ingenuo sobre el mismo tramo y el mismo h
  ing <- ajustar_ingenuo(y_est, s = s)
  pronos_ing <- ing$pronosticar(hh)
  esc <- escala_mase(y_est, s = s)
  
  med_dentro <- medidas(y_est, ajuste$yhat, escala = esc)
  med_fuera  <- medidas(y_val, pronos, escala = esc)
  med_ing    <- medidas(y_val, pronos_ing, escala = esc)
  
  # 4. Validacion de errores de un paso
  val <- validar_errores(e_est, fechas = part$estimacion$fecha, p = p,
                         titulo = sprintf("%s con %s: errores de un paso", nombre, etiqueta),
                         imprimir = imprimir)
  
  # 5. Grafico final
  g_final <- graficar_ajuste(datos, ajuste$yhat, pronos, part$n_estimacion,
                             sprintf("%s: ajuste y %d pronosticos con %s", nombre, hh, etiqueta),
                             attr(datos, "unidad"), attr(datos, "fuente"))
  
  if (imprimir) {
    cat("\n########## Ejemplo:", nombre, "con", etiqueta, "##########\n")
    print(as.data.frame(descripcion))
    cat(sprintf("\nParticion: estimacion n = %d, validacion h = %d.\n", part$n_estimacion, hh))
    cat("\nLjung-Box sobre la serie (p = 0):\n")
    reportar_prueba(lb_serie)
    cat("\nRezagos de la serie que exceden la banda de ruido blanco: ",
        paste(rh_serie$rezago[rh_serie$significativo], collapse = ", "), "\n", sep = "")
    cat("\nMedidas de un paso dentro del tramo de estimacion:\n")
    print(as.data.frame(med_dentro))
    cat("\nMedidas de los", hh, "pronosticos sobre el tramo de validacion:\n")
    print(as.data.frame(med_fuera))
    cat("\nReferente", ing$parametros$metodo, "sobre el mismo tramo:\n")
    print(as.data.frame(med_ing))
  }
  
  list(
    nombre = nombre, etiqueta = etiqueta, h = hh, p = p,
    particion = part, descripcion = descripcion,
    grafico_serie = g_serie, correlograma_serie = corr_serie,
    rh_serie = rh_serie, ljung_box_serie = lb_serie,
    ajuste = ajuste, errores = e_est, pronosticos = pronos,
    referente = ing, pronosticos_referente = pronos_ing, escala_mase = esc,
    medidas_dentro = med_dentro, medidas_fuera = med_fuera, medidas_referente = med_ing,
    validacion = val, grafico_final = g_final
  )
}

#' Guarda una figura en figs/ con dimensiones uniformes.
guardar_fig <- function(grafico, archivo, ancho = 9, alto = 5.5) {
  ruta <- file.path("figs", archivo)
  ggplot2::ggsave(ruta, grafico, width = ancho, height = alto, dpi = 150)
  invisible(ruta)
}