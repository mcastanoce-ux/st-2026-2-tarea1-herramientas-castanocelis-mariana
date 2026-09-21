suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

# Convierte el tiempo decimal de un objeto ts a una fecha de clase Date.
.tiempo_a_fecha <- function(tt, frecuencia) {
  stopifnot(is.numeric(tt), length(frecuencia) == 1L)
  anio <- floor(tt + 1e-8)
  if (frecuencia == 12) {
    mes <- round((tt - anio) * 12) + 1
    as.Date(sprintf("%04d-%02d-01", anio, mes))
  } else if (frecuencia == 4) {
    trimestre <- round((tt - anio) * 4) + 1
    as.Date(sprintf("%04d-%02d-01", anio, (trimestre - 1) * 3 + 1))
  } else if (frecuencia == 1) {
    as.Date(sprintf("%04d-01-01", anio))
  } else {
    stop("leer_serie(): solo se admiten frecuencias 1, 4 o 12.", call. = FALSE)
  }
}

# Sirve para verificar equiespaciamiento.
.indice_periodo <- function(fecha, frecuencia) {
  anio <- as.integer(format(fecha, "%Y"))
  mes  <- as.integer(format(fecha, "%m"))
  if (frecuencia == 12) {
    anio * 12L + (mes - 1L)
  } else if (frecuencia == 4) {
    anio * 4L + ((mes - 1L) %/% 3L)
  } else {
    anio
  }
}


leer_serie <- function(x, fuente, unidad, frecuencia = NULL) {
  stopifnot(
    "fuente debe ser una cadena de texto" = is.character(fuente) && length(fuente) == 1L,
    "unidad debe ser una cadena de texto" = is.character(unidad) && length(unidad) == 1L
  )
  
  if (stats::is.ts(x)) {
    stopifnot("la serie ts debe ser univariada" = is.null(dim(x)) || ncol(x) == 1L)
    frecuencia <- stats::frequency(x)
    valores <- as.numeric(x)
    fechas  <- .tiempo_a_fecha(as.numeric(stats::time(x)), frecuencia)
  } else if (is.character(x) && length(x) == 1L) {
    stopifnot("el archivo csv no existe" = file.exists(x))
    bruto <- utils::read.csv(x, stringsAsFactors = FALSE)
    stopifnot(
      "el csv debe tener columnas 'fecha' y 'valor'" =
        all(c("fecha", "valor") %in% names(bruto))
    )
    fechas  <- as.Date(bruto$fecha)
    valores <- as.numeric(bruto$valor)
    stopifnot("hay fechas que no se pudieron convertir a Date" = !any(is.na(fechas)))
    if (is.null(frecuencia)) {
      paso <- stats::median(as.numeric(diff(fechas)))
      frecuencia <- if (paso > 300) 1 else if (paso > 75) 4 else 12
    }
  } else {
    stop("leer_serie(): x debe ser un objeto ts o la ruta a un archivo csv.",
         call. = FALSE)
  }
  
  stopifnot(
    "la serie no puede contener NA" = !any(is.na(valores)),
    "la serie debe tener al menos 3 observaciones" = length(valores) >= 3L,
    "la frecuencia debe ser 1, 4 o 12" = frecuencia %in% c(1, 4, 12)
  )
  
  # Fechas crecientes
  if (any(diff(fechas) <= 0)) {
    stop("leer_serie(): las fechas no son estrictamente crecientes.", call. = FALSE)
  }
  # Fechas equiespaciadas segun la frecuencia declarada
  saltos <- diff(.indice_periodo(fechas, frecuencia))
  if (any(saltos != 1L)) {
    malas <- which(saltos != 1L)[1]
    stop(sprintf(
      "leer_serie(): las fechas no son equiespaciadas para frecuencia %g. Primer hueco entre %s y %s.",
      frecuencia, fechas[malas], fechas[malas + 1L]), call. = FALSE)
  }
  
  datos <- tibble::tibble(
    t     = seq_along(valores),
    fecha = fechas,
    y     = valores
  )
  attr(datos, "frecuencia") <- frecuencia
  attr(datos, "fuente")     <- fuente
  attr(datos, "unidad")     <- unidad
  datos
}


#' rango de fechas y numero de observaciones.
describir_serie <- function(datos, nombre) {
  stopifnot(
    "datos debe venir de leer_serie()" = all(c("t", "fecha", "y") %in% names(datos)),
    "nombre debe ser una cadena" = is.character(nombre) && length(nombre) == 1L
  )
  tibble::tibble(
    serie          = nombre,
    fuente         = attr(datos, "fuente"),
    unidad         = attr(datos, "unidad"),
    frecuencia     = attr(datos, "frecuencia"),
    inicio         = format(min(datos$fecha)),
    fin            = format(max(datos$fecha)),
    observaciones  = nrow(datos),
    media          = mean(datos$y),
    desv_estandar  = stats::sd(datos$y),
    minimo         = min(datos$y),
    maximo         = max(datos$y)
  )
}

#' Particion en tramo de estimacion y tramo de validacion.

particionar <- function(datos, h = NULL, estacional = FALSE) {
  stopifnot("datos debe venir de leer_serie()" = all(c("t", "fecha", "y") %in% names(datos)))
  T_total <- nrow(datos)
  s <- attr(datos, "frecuencia")
  if (is.null(h)) {
    h <- min(12, floor(0.2 * T_total))
    if (estacional && s > 1) h <- max(h, s)
  }
  stopifnot(
    "h debe ser un entero positivo menor que T" = h >= 1 && h < T_total,
    "el tramo de estimacion queda demasiado corto" = (T_total - h) >= 10
  )
  n_est <- T_total - h
  # El subconjunto de un tibble pierde los atributos propios; se vuelven a pegar.
  pegar_attr <- function(d) {
    attr(d, "frecuencia") <- attr(datos, "frecuencia")
    attr(d, "fuente")     <- attr(datos, "fuente")
    attr(d, "unidad")     <- attr(datos, "unidad")
    d
  }
  list(
    estimacion = pegar_attr(datos[seq_len(n_est), ]),
    validacion = pegar_attr(datos[(n_est + 1L):T_total, ]),
    h = h,
    n_estimacion = n_est
  )
}