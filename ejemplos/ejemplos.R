set.seed(2026)   # se fija por reproducibilidad

for (archivo in c("R/00-lectura.R", "R/01-graficos.R",
                  "R/02-metodos.R", "R/03-evaluacion.R")) {
  if (!file.exists(archivo)) {
    stop("ejemplos.R debe ejecutarse desde la raiz del repositorio. No se encuentra ",
         archivo, call. = FALSE)
  }
  source(archivo)
}
if (!dir.exists("figs")) dir.create("figs")

resumen_global <- list()


# BLOQUE DE VERIFICACION

cat("\n==========================================================\n")
cat(" BLOQUE DE VERIFICACION\n")
cat("==========================================================\n")
for (nombre_serie in c("Nile", "LakeHuron", "austres", "WWWusage")) {
  yv <- as.numeric(get(nombre_serie))
  mv <- min(floor(length(yv) / 4), 24)
  r_propio <- acf_manual(yv, mv)$r
  r_R <- as.numeric(stats::acf(yv, lag.max = mv, plot = FALSE)$acf)[-1]
  dif_acf <- max(abs(r_propio - r_R))
  q_propio <- ljung_box(r_propio, T = length(yv), m = mv, p = 0)
  q_R <- stats::Box.test(yv, lag = mv, type = "Ljung-Box")
  cat(sprintf("%-12s  max|r_h propio - acf()| = %.3e (< 1e-12: %s) | Q_m propio = %.6f vs Box.test = %.6f, dif = %.3e\n",
              nombre_serie, dif_acf, dif_acf < 1e-12,
              q_propio$estadistico, as.numeric(q_R$statistic),
              abs(q_propio$estadistico - as.numeric(q_R$statistic))))
}
# Coeficientes de las tendencias contra lm()
y_ver <- as.numeric(austres); t_ver <- seq_along(y_ver)
dif_lineal <- max(abs(ajustar_tendencia(y_ver, "lineal")$parametros$coeficientes -
                        coef(stats::lm(y_ver ~ t_ver))))
dif_cuad <- max(abs(ajustar_tendencia(y_ver, "cuadratica")$parametros$coeficientes -
                      coef(stats::lm(y_ver ~ t_ver + I(t_ver^2)))))
y_exp <- as.numeric(JohnsonJohnson); t_exp <- seq_along(y_exp)
dif_exp <- max(abs(ajustar_tendencia(y_exp, "exponencial")$parametros$coeficientes -
                     coef(stats::lm(log(y_exp) ~ t_exp))))
cat(sprintf("Coeficientes propios vs lm(): lineal %.3e | cuadratica %.3e | exponencial (en logs) %.3e\n",
            dif_lineal, dif_cuad, dif_exp))
# SES: forma de correccion de error vs promedio ponderado; Holt: correccion de error
cat(sprintf("SES: max|correccion de error - promedio ponderado| = %.3e\n",
            ajustar_ses(as.numeric(Nile), 0.3)$parametros$dif_max_formas))
cat(sprintf("Holt: max|ecuaciones - forma de correccion de error| = %.3e\n",
            ajustar_holt(as.numeric(WWWusage), 0.5, 0.3)$parametros$dif_max_correccion_error))

cat("\n==========================================================\n")
cat(" EJEMPLO 1. Media simple sobre Nile\n")
cat("==========================================================\n")
# Nile: caudal anual del Nilo en Ashwan. Serie de nivel sin tendencia sostenida;
# el supuesto de la media simple (nivel constante, sin tendencia) es defendible.
nile <- leer_serie(Nile,
                   fuente = "R datasets: Nile. Durbin y Koopman (2001), medicion anual del caudal del rio Nilo en Ashwan, 1871-1970",
                   unidad = "Caudal anual (10^8 m^3)")
ej1 <- protocolo_ejemplo(nile, "Nile", ajustar_media, "media simple", p = 1)
guardar_fig(ej1$grafico_serie,        "ej1-nile-serie.png")
guardar_fig(ej1$correlograma_serie,   "ej1-nile-correlograma.png", alto = 6.5)
guardar_fig(ej1$validacion$grafico_errores, "ej1-nile-errores.png", alto = 4)
guardar_fig(ej1$validacion$correlograma,    "ej1-nile-correlograma-errores.png", alto = 6.5)
guardar_fig(ej1$grafico_final,        "ej1-nile-ajuste.png")
resumen_global$ej1 <- ej1

cat("\n==========================================================\n")
cat(" EJEMPLO 2. Media movil sobre LakeHuron\n")
cat("==========================================================\n")
# LakeHuron: nivel anual del lago Huron. Nivel que deriva lentamente, sin
# tendencia lineal marcada: es el caso de la media movil.
huron <- leer_serie(LakeHuron,
                    fuente = "R datasets: LakeHuron. Brockwell y Davis (1991), nivel anual del lago Huron en julio, 1875-1972",
                    unidad = "Nivel del lago (pies)")
part2 <- particionar(huron)
opt2 <- optimizar(part2$estimacion$y, "mm")
cat("\nOptimizacion de k (media movil):\n"); print(as.data.frame(opt2$rejilla))
cat("Optimo: k =", opt2$optimo$k, " MSE =", opt2$optimo$mse,
    if (opt2$en_borde) "(en el borde de la rejilla)" else "", "\n")
guardar_fig(graficar_optimizacion(opt2, "LakeHuron: MSE de un paso contra k"),
            "ej2-huron-optimizacion.png", alto = 4)
k2 <- opt2$optimo$k
ej2 <- protocolo_ejemplo(huron, "LakeHuron", function(y) ajustar_mm(y, k2),
                         sprintf("media movil (k = %d)", k2), p = 1)
guardar_fig(ej2$grafico_serie,      "ej2-huron-serie.png")
guardar_fig(ej2$correlograma_serie, "ej2-huron-correlograma.png", alto = 6.5)
guardar_fig(ej2$validacion$grafico_errores, "ej2-huron-errores.png", alto = 4)
guardar_fig(ej2$validacion$correlograma,    "ej2-huron-correlograma-errores.png", alto = 6.5)
guardar_fig(ej2$grafico_final,      "ej2-huron-ajuste.png")
resumen_global$ej2 <- ej2

cat("\n==========================================================\n")
cat(" EJEMPLO 3. Suavizamiento exponencial simple sobre Nile\n")
cat("==========================================================\n")
part3 <- particionar(nile)
opt3 <- optimizar(part3$estimacion$y, "ses")
cat("\nOptimo alpha =", opt3$optimo$alpha, " MSE =", opt3$optimo$mse,
    if (opt3$en_borde) "(en el borde de la rejilla)" else "", "\n")
guardar_fig(graficar_optimizacion(opt3, "Nile: MSE de un paso contra alpha (SES)"),
            "ej3-nile-optimizacion.png", alto = 4)
a3 <- opt3$optimo$alpha
ej3 <- protocolo_ejemplo(nile, "Nile", function(y) ajustar_ses(y, a3),
                         sprintf("SES (alpha = %.2f)", a3), p = 1)
guardar_fig(ej3$validacion$grafico_errores, "ej3-nile-ses-errores.png", alto = 4)
guardar_fig(ej3$validacion$correlograma,    "ej3-nile-ses-correlograma-errores.png", alto = 6.5)
guardar_fig(ej3$grafico_final,              "ej3-nile-ses-ajuste.png")
resumen_global$ej3 <- ej3

cat("\n==========================================================\n")
cat(" EJEMPLO 4. Doble media movil sobre austres\n")
cat("==========================================================\n")
# austres: residentes en Australia, trimestral. Tendencia creciente muy regular
# y sin estacionalidad apreciable: es el caso de la doble media movil.
austres_d <- leer_serie(austres,
                        fuente = "R datasets: austres. Brockwell y Davis (1996), numero de residentes en Australia, trimestral 1971:2-1993:3",
                        unidad = "Residentes (miles de personas)")
part4 <- particionar(austres_d)
opt4 <- optimizar(part4$estimacion$y, "dmm")
cat("\nOptimo k (doble media movil) =", opt4$optimo$k, " MSE =", opt4$optimo$mse,
    if (opt4$en_borde) "(en el borde de la rejilla)" else "", "\n")
guardar_fig(graficar_optimizacion(opt4, "austres: MSE de un paso contra k (doble media movil)"),
            "ej4-austres-optimizacion.png", alto = 4)
k4 <- opt4$optimo$k
ej4 <- protocolo_ejemplo(austres_d, "austres", function(y) ajustar_dmm(y, k4),
                         sprintf("doble media movil (k = %d)", k4), p = 2)
guardar_fig(ej4$grafico_serie,      "ej4-austres-serie.png")
guardar_fig(ej4$correlograma_serie, "ej4-austres-correlograma.png", alto = 6.5)
guardar_fig(ej4$validacion$grafico_errores, "ej4-austres-errores.png", alto = 4)
guardar_fig(ej4$validacion$correlograma,    "ej4-austres-correlograma-errores.png", alto = 6.5)
guardar_fig(ej4$grafico_final,      "ej4-austres-ajuste.png")
resumen_global$ej4 <- ej4

cat("\n==========================================================\n")
cat(" EJEMPLO 5. Tendencia lineal sobre austres\n")
cat("==========================================================\n")
ej5 <- protocolo_ejemplo(austres_d, "austres",
                         function(y) ajustar_tendencia(y, "lineal"),
                         "tendencia lineal", p = 2)
cat("\nTabla de coeficientes (tendencia lineal, austres):\n")
print(as.data.frame(ej5$ajuste$parametros$tabla_coeficientes))
cat(sprintf("R2 = %.4f ; sigma2 = %.4f ; DW = %.4f ; rezagos HAC = %d ; gl = %d\n",
            ej5$ajuste$parametros$R2, ej5$ajuste$parametros$sigma2,
            ej5$ajuste$parametros$durbin_watson,
            ej5$ajuste$parametros$rezagos_hac, ej5$ajuste$parametros$gl))
for (j in seq_len(nrow(ej5$ajuste$parametros$tabla_coeficientes))) {
  tc <- ej5$ajuste$parametros$tabla_coeficientes[j, ]
  reportar_prueba(t_coeficiente(tc$estimacion, tc$ee_hac,
                                ej5$ajuste$parametros$gl, tc$coeficiente))
}
guardar_fig(ej5$validacion$grafico_errores, "ej5-austres-lineal-errores.png", alto = 4)
guardar_fig(ej5$validacion$correlograma,    "ej5-austres-lineal-correlograma-errores.png", alto = 6.5)
guardar_fig(ej5$grafico_final,              "ej5-austres-lineal-ajuste.png")
# Residuos contra valores ajustados (varianza constante)
g5_res <- ggplot(tibble::tibble(ajustado = ej5$ajuste$yhat, residuo = ej5$errores),
                 aes(x = ajustado, y = residuo)) +
  geom_hline(yintercept = 0, colour = "grey40") + geom_point(size = 1) +
  labs(title = "austres, tendencia lineal: residuos contra valores ajustados",
       x = "Valor ajustado", y = "Residuo") + theme_bw(base_size = 10)
guardar_fig(g5_res, "ej5-austres-lineal-residuos-ajustados.png", alto = 4)
resumen_global$ej5 <- ej5

cat("\n==========================================================\n")
cat(" EJEMPLO 6. Tendencia cuadratica sobre airmiles\n")
cat("==========================================================\n")
# airmiles: millas-pasajero de las aerolineas comerciales de EE. UU. El
# crecimiento es convexo: la tendencia cuadratica es el candidato natural.
airm <- leer_serie(airmiles,
                   fuente = "R datasets: airmiles. Brown (1963), millas-pasajero de las aerolineas comerciales de EE. UU., 1937-1960",
                   unidad = "Millas-pasajero (miles)")
ej6 <- protocolo_ejemplo(airm, "airmiles",
                         function(y) ajustar_tendencia(y, "cuadratica"),
                         "tendencia cuadratica", p = 3)
cat("\nTabla de coeficientes (tendencia cuadratica, airmiles):\n")
print(as.data.frame(ej6$ajuste$parametros$tabla_coeficientes))
cat(sprintf("R2 = %.4f ; sigma2 = %.4f ; DW = %.4f ; rezagos HAC = %d ; gl = %d\n",
            ej6$ajuste$parametros$R2, ej6$ajuste$parametros$sigma2,
            ej6$ajuste$parametros$durbin_watson,
            ej6$ajuste$parametros$rezagos_hac, ej6$ajuste$parametros$gl))
for (j in seq_len(nrow(ej6$ajuste$parametros$tabla_coeficientes))) {
  tc <- ej6$ajuste$parametros$tabla_coeficientes[j, ]
  reportar_prueba(t_coeficiente(tc$estimacion, tc$ee_hac,
                                ej6$ajuste$parametros$gl, tc$coeficiente))
}
guardar_fig(ej6$grafico_serie,      "ej6-airmiles-serie.png")
guardar_fig(ej6$correlograma_serie, "ej6-airmiles-correlograma.png", alto = 6.5)
guardar_fig(ej6$validacion$grafico_errores, "ej6-airmiles-errores.png", alto = 4)
guardar_fig(ej6$validacion$correlograma,    "ej6-airmiles-correlograma-errores.png", alto = 6.5)
guardar_fig(ej6$grafico_final,      "ej6-airmiles-ajuste.png")
resumen_global$ej6 <- ej6

cat("\n==========================================================\n")
cat(" EJEMPLO 7. Tendencia exponencial sobre JohnsonJohnson\n")
cat("==========================================================\n")
# JohnsonJohnson: utilidades trimestrales por accion. Crecimiento multiplicativo
# claro; se declara como limite la estacionalidad trimestral, que el metodo no
# modela y que reaparece en el correlograma de los errores.
jj <- leer_serie(JohnsonJohnson,
                 fuente = "R datasets: JohnsonJohnson. Shumway y Stoffer (2000), utilidades trimestrales por accion de Johnson & Johnson, 1960-1980",
                 unidad = "Utilidad por accion (dolares)")
ej7 <- protocolo_ejemplo(jj, "JohnsonJohnson",
                         function(y) ajustar_tendencia(y, "exponencial", corregir_sesgo = TRUE),
                         "tendencia exponencial", p = 2, estacional = TRUE)
cat("\nTabla de coeficientes en la escala logaritmica (tendencia exponencial, JohnsonJohnson):\n")
print(as.data.frame(ej7$ajuste$parametros$tabla_coeficientes))
cat(sprintf("b0 = exp(a) = %.6f ; b1 = exp(theta) = %.6f ; factor de correccion de sesgo = %.6f\n",
            ej7$ajuste$parametros$b0_exponencial, ej7$ajuste$parametros$b1_exponencial,
            ej7$ajuste$parametros$factor_sesgo))
cat(sprintf("R2 (en logaritmos) = %.4f ; sigma2_ln = %.4f ; DW = %.4f ; gl = %d\n",
            ej7$ajuste$parametros$R2, ej7$ajuste$parametros$sigma2,
            ej7$ajuste$parametros$durbin_watson, ej7$ajuste$parametros$gl))
for (j in seq_len(nrow(ej7$ajuste$parametros$tabla_coeficientes))) {
  tc <- ej7$ajuste$parametros$tabla_coeficientes[j, ]
  reportar_prueba(t_coeficiente(tc$estimacion, tc$ee_hac,
                                ej7$ajuste$parametros$gl, tc$coeficiente))
}
guardar_fig(ej7$grafico_serie,      "ej7-jj-serie.png")
guardar_fig(ej7$correlograma_serie, "ej7-jj-correlograma.png", alto = 6.5)
guardar_fig(ej7$validacion$grafico_errores, "ej7-jj-errores.png", alto = 4)
guardar_fig(ej7$validacion$correlograma,    "ej7-jj-correlograma-errores.png", alto = 6.5)
guardar_fig(ej7$grafico_final,      "ej7-jj-ajuste.png")
resumen_global$ej7 <- ej7

cat("\n==========================================================\n")
cat(" EJEMPLO 8. Holt lineal sobre WWWusage\n")
cat("==========================================================\n")
# WWWusage: usuarios conectados a un servidor, minuto a minuto. Nivel y
# pendiente locales que cambian en el tiempo, sin estacionalidad: caso de Holt.
www <- leer_serie(WWWusage,
                  fuente = "R datasets: WWWusage. Durbin y Koopman (2001), numero de usuarios conectados a un servidor de internet, un dato por minuto",
                  unidad = "Usuarios conectados")
part8 <- particionar(www)
opt8 <- optimizar(part8$estimacion$y, "holt")
cat("\nOptimo Holt: alpha =", opt8$optimo$alpha, " beta =", opt8$optimo$beta,
    " MSE =", opt8$optimo$mse,
    if (opt8$en_borde) "(en el borde de la rejilla)" else "", "\n")
guardar_fig(graficar_optimizacion(opt8, "WWWusage: mapa de MSE sobre la rejilla (alpha, beta)"),
            "ej8-www-optimizacion.png", alto = 5)
a8 <- opt8$optimo$alpha; b8 <- opt8$optimo$beta
ej8 <- protocolo_ejemplo(www, "WWWusage", function(y) ajustar_holt(y, a8, b8),
                         sprintf("Holt lineal (alpha = %.2f, beta = %.2f)", a8, b8), p = 2)
guardar_fig(ej8$grafico_serie,      "ej8-www-serie.png")
guardar_fig(ej8$correlograma_serie, "ej8-www-correlograma.png", alto = 6.5)
guardar_fig(ej8$validacion$grafico_errores, "ej8-www-errores.png", alto = 4)
guardar_fig(ej8$validacion$correlograma,    "ej8-www-correlograma-errores.png", alto = 6.5)
guardar_fig(ej8$grafico_final,      "ej8-www-ajuste.png")
resumen_global$ej8 <- ej8

cat("\n==========================================================\n")
cat(" CONTRAEJEMPLO (4e). Media simple sobre AirPassengers\n")
cat("==========================================================\n")
# AirPassengers tiene tendencia creciente y estacionalidad multiplicativa: el
# supuesto de nivel constante de la media simple se viola por completo.
aire <- leer_serie(AirPassengers,
                   fuente = "R datasets: AirPassengers. Box, Jenkins y Reinsel (1976), pasajeros mensuales de aerolineas internacionales, 1949-1960",
                   unidad = "Pasajeros (miles)")
ejc <- protocolo_ejemplo(aire, "AirPassengers", ajustar_media,
                         "media simple (contraejemplo)", p = 1, estacional = TRUE)
guardar_fig(ejc$grafico_serie,      "ejc-airpassengers-serie.png")
guardar_fig(ejc$correlograma_serie, "ejc-airpassengers-correlograma.png", alto = 6.5)
guardar_fig(ejc$validacion$grafico_errores, "ejc-airpassengers-errores.png", alto = 4)
guardar_fig(ejc$validacion$correlograma,    "ejc-airpassengers-correlograma-errores.png", alto = 6.5)
guardar_fig(ejc$grafico_final,      "ejc-airpassengers-ajuste.png")
resumen_global$ejc <- ejc

# ---------------------------------------------------------------------------
# Tabla resumen de resultados (una fila por ejemplo) para el README
# ---------------------------------------------------------------------------
fila <- function(ej, parametros) {
  tibble::tibble(
    ejemplo = ej$etiqueta,
    serie = ej$nombre,
    parametros = parametros,
    h = ej$h,
    MASE_metodo = ej$medidas_fuera$MASE,
    MASE_referente = ej$medidas_referente$MASE,
    MSE_validacion = ej$medidas_fuera$MSE,
    MSE_referente = ej$medidas_referente$MSE
  )
}
tabla_resumen <- dplyr::bind_rows(
  fila(ej1, "ninguno"),
  fila(ej2, sprintf("k = %d", k2)),
  fila(ej3, sprintf("alpha = %.2f", a3)),
  fila(ej4, sprintf("k = %d", k4)),
  fila(ej5, "MCO, 2 coeficientes"),
  fila(ej6, "MCO, 3 coeficientes"),
  fila(ej7, "MCO en logaritmos, correccion de sesgo"),
  fila(ej8, sprintf("alpha = %.2f, beta = %.2f", a8, b8)),
  fila(ejc, "contraejemplo")
)
cat("\n==========================================================\n")
cat(" RESUMEN DE RESULTADOS\n")
cat("==========================================================\n")
print(as.data.frame(tabla_resumen), digits = 4)
utils::write.csv(tabla_resumen, "figs/tabla-resumen.csv", row.names = FALSE)

cat("\nEjecucion terminada. Figuras guardadas en figs/.\n")
