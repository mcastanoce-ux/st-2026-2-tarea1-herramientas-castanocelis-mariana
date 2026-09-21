
suppressPackageStartupMessages({
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(patchwork)
})


acf_manual <- function(y, m = NULL) {
  stopifnot(
    "y debe ser numerico" = is.numeric(y),
    "y no puede contener NA" = !any(is.na(y)),
    "y debe tener al menos 4 observaciones" = length(y) >= 4L
  )
  n <- length(y)
  if (is.null(m)) m <- min(floor(n / 4), 24)
  stopifnot(
    "m debe ser un entero positivo menor que n" = m >= 1 && m < n
  )
  d <- y - mean(y)
  c0 <- sum(d^2) / n                      # autocovarianza de rezago 0
  r <- vapply(seq_len(m), function(h) {
    (sum(d[(h + 1):n] * d[1:(n - h)]) / n) / c0
  }, numeric(1))
  tibble::tibble(rezago = seq_len(m), r = r)
}
#PACF
pacf_stats <- function(y, m = NULL) {
  stopifnot("y debe ser numerico sin NA" = is.numeric(y) && !any(is.na(y)))
  n <- length(y)
  if (is.null(m)) m <- min(floor(n / 4), 24)
  sal <- stats::pacf(y, lag.max = m, plot = FALSE)
  tibble::tibble(rezago = seq_len(m), phi = as.numeric(sal$acf)[seq_len(m)])
}

#Banda de ruido blanco 
banda_ruido_blanco <- function(n, ci = 0.95) {
  stopifnot("n debe ser un entero positivo" = is.numeric(n) && n > 0)
  stats::qnorm((1 + ci) / 2) / sqrt(n)
}

#Grafico de la serie en el tiempo.
graficar_serie <- function(datos, titulo) {
  stopifnot(
    "datos debe venir de leer_serie()" = all(c("t", "fecha", "y") %in% names(datos)),
    "titulo debe ser una cadena" = is.character(titulo) && length(titulo) == 1L
  )
  unidad <- attr(datos, "unidad")
  fuente <- attr(datos, "fuente")
  frec   <- attr(datos, "frecuencia")
  etiqueta_frec <- c("1" = "anual", "4" = "trimestral", "12" = "mensual")[as.character(frec)]
  
  ggplot(datos, aes(x = fecha, y = y)) +
    geom_line(linewidth = 0.6, colour = "grey20") +
    geom_point(size = 0.7, colour = "grey20") +
    scale_x_date(date_labels = "%Y", guide = guide_axis(check.overlap = TRUE)) +
    labs(
      title    = titulo,
      subtitle = sprintf("Frecuencia %s (%g observaciones por anio)", etiqueta_frec, frec),
      x        = "Fecha",
      y        = unidad,
      caption  = sprintf("Fuente: %s. n = %d observaciones (%s a %s).",
                         fuente, nrow(datos),
                         format(min(datos$fecha)), format(max(datos$fecha)))
    ) +
    theme_bw(base_size = 11) +
    theme(plot.caption = element_text(hjust = 0, size = 8))
}


correlograma <- function(datos, m = NULL, titulo = "Correlograma", verificar = TRUE) {
  y <- if (is.numeric(datos)) datos else datos$y
  stopifnot(
    "la sucesion no puede contener NA" = !any(is.na(y)),
    "la sucesion debe tener al menos 4 observaciones" = length(y) >= 4L
  )
  n <- length(y)
  if (is.null(m)) m <- min(floor(n / 4), 24)
  stopifnot("m debe ser un entero positivo menor que n" = m >= 1 && m < n)
  
  tabla_acf  <- acf_manual(y, m)
  tabla_pacf <- pacf_stats(y, m)
  banda <- banda_ruido_blanco(n)          # n = longitud de la sucesion graficada
  
  # Bloque de verificacion: ACF propia contra stats::acf()
  verif <- NA_real_
  if (verificar) {
    ref <- as.numeric(stats::acf(y, lag.max = m, plot = FALSE)$acf)[-1]
    verif <- max(abs(tabla_acf$r - ref))
  }
  
  base_tema <- theme_bw(base_size = 10) +
    theme(plot.caption = element_text(hjust = 0, size = 8))
  
  g_acf <- ggplot(tabla_acf, aes(x = rezago, y = r)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-banda, banda), linetype = "dashed", colour = "blue") +
    geom_segment(aes(xend = rezago, yend = 0), linewidth = 0.7) +
    scale_x_continuous(breaks = scales_enteros(m)) +
    labs(title = titulo,
         subtitle = sprintf("ACF muestral con divisor T; banda +-1.96/sqrt(n), n = %d", n),
         x = "Rezago h", y = expression(r[h])) +
    base_tema
  
  g_pacf <- ggplot(tabla_pacf, aes(x = rezago, y = phi)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-banda, banda), linetype = "dashed", colour = "blue") +
    geom_segment(aes(xend = rezago, yend = 0), linewidth = 0.7) +
    scale_x_continuous(breaks = scales_enteros(m)) +
    labs(subtitle = "PACF muestral (stats::pacf, plot = FALSE)",
         x = "Rezago h", y = expression(hat(phi)[hh]),
         caption = if (verificar)
           sprintf("Verificacion: max |r_h propio - acf()| = %.3e (debe ser < 1e-12).", verif)
         else NULL) +
    base_tema
  
  panel <- g_acf / g_pacf
  attr(panel, "tabla") <- tabla_acf
  attr(panel, "pacf") <- tabla_pacf
  attr(panel, "banda") <- banda
  attr(panel, "verificacion") <- verif
  panel
}

# Cortes enteros y legibles para el eje de rezagos.
scales_enteros <- function(m) {
  paso <- if (m <= 12) 1 else if (m <= 24) 2 else 4
  seq(paso, m, by = paso)
}

#' Contraste individual de los r_h contra la banda de ruido blanco.
contrastar_rh <- function(y, m = NULL) {
  y <- if (is.numeric(y)) y else y$y
  n <- length(y)
  if (is.null(m)) m <- min(floor(n / 4), 24)
  tabla <- acf_manual(y, m)
  banda <- banda_ruido_blanco(n)
  tabla %>%
    dplyr::mutate(
      banda = banda,
      z = r * sqrt(n),                        # estadistico z_h = r_h * sqrt(n)
      significativo = abs(r) > banda,
      valor_p = 2 * (1 - stats::pnorm(abs(z)))
    )
}

#' Grafico final de un ejemplo: serie observada, ajuste dentro del tramo de
graficar_ajuste <- function(datos, yhat_est, pronos, corte, titulo, unidad, fuente) {
  stopifnot("corte debe ser menor que el numero de observaciones" = corte < nrow(datos))
  h <- length(pronos)
  serie <- dplyr::mutate(datos, tramo = ifelse(t <= corte, "estimacion", "validacion"))
  ajuste <- tibble::tibble(fecha = datos$fecha[seq_len(corte)], y = yhat_est)
  futuro <- tibble::tibble(fecha = datos$fecha[(corte + 1):(corte + h)], y = pronos)
  
  ggplot(serie, aes(x = fecha, y = y)) +
    geom_line(colour = "grey25", linewidth = 0.6) +
    geom_line(data = ajuste, aes(x = fecha, y = y), colour = "#1b6ca8",
              linewidth = 0.6, na.rm = TRUE) +
    geom_line(data = futuro, aes(x = fecha, y = y), colour = "#c0392b", linewidth = 0.8) +
    geom_point(data = futuro, aes(x = fecha, y = y), colour = "#c0392b", size = 1.2) +
    geom_vline(xintercept = as.numeric(datos$fecha[corte]), linetype = "dotted") +
    labs(title = titulo,
         subtitle = "Gris: serie observada. Azul: ajuste de un paso. Rojo: pronosticos extramuestrales.",
         x = "Fecha", y = unidad,
         caption = sprintf("Fuente: %s. Linea punteada: fin del tramo de estimacion (h = %d).",
                           fuente, h)) +
    theme_bw(base_size = 11) +
    theme(plot.caption = element_text(hjust = 0, size = 8))
}
