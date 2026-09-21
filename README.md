https://github.com/mcastanoce-ux/st-2026-2-tarea1-herramientas-castanocelis-mariana

# Tarea 1: caja de herramientas de pronóstico

Series de tiempo univariadas , semestre 2026-II. Universidad Nacional
de Colombia, sede Medellín. Prof. Juan Pablo Valencia Arango.
Autor: Mariana Castaño Celis.

## 1. Qué contiene el repositorio

| Archivo | Contenido | Dependencias |
|---|---|---|
| `R/00-lectura.R` | `leer_serie()`, `describir_serie()`, `particionar()` | tibble, dplyr |
| `R/01-graficos.R` | `graficar_serie()`, `correlograma()`, `acf_manual()`, `pacf_stats()`, `banda_ruido_blanco()`, `contrastar_rh()`, `graficar_ajuste()` | ggplot2, patchwork, tibble, dplyr, stats |
| `R/02-metodos.R` | `ajustar_media()`, `ajustar_mm()`, `ajustar_ses()`, `ajustar_dmm()`, `ajustar_tendencia()`, `ajustar_holt()`, `optimizar()`, `ajustar_ingenuo()` | tibble, dplyr, R base |
| `R/03-evaluacion.R` | `medidas()`, `escala_mase()`, `ljung_box()`, `jarque_bera()`, `durbin_watson()`, `t_media_cero()`, `t_coeficiente()`, `reportar_prueba()`, `validar_errores()`, `graficar_optimizacion()`, `protocolo_ejemplo()`, `guardar_fig()` | ggplot2, patchwork, tibble, dplyr, stats |
| `ejemplos/ejemplos.R` | Carga `R/` y ejecuta los ocho ejemplos más el contraejemplo | las anteriores |
| `informe/informe.qmd` | Análisis detallado de los ocho ejemplos | quarto, knitr |
| `informe/informe.html` | El `.qmd` renderizado | — |
| `figs/` | Figuras generadas por `ejemplos.R` (se versionan) | — |
| `sesion-info.txt` | Salida de `sessionInfo()` al cerrar la entrega | — |

Ningún archivo usa `forecast`, `fable`, `smooth`, `TTR` ni `zoo`. De `stats` no
se usan `HoltWinters()`, `filter()`, `decompose()` ni `stl()`; `acf()`, `pacf()`,
`lm()` y `Box.test()` aparecen únicamente en bloques de verificación, y `pacf()`
además para dibujar la PACF, como autoriza el enunciado.

## 2. Cómo se corre

Desde la raíz del repositorio, con R 4.3 o posterior:

```r
install.packages(c("dplyr", "tidyr", "purrr", "tibble", "ggplot2", "patchwork"))
source("ejemplos/ejemplos.R")
```

Para el informe:

```bash
quarto render informe/informe.qmd
```



## 3. Cómo se usan las funciones

```r
source("R/00-lectura.R"); source("R/01-graficos.R")
source("R/02-metodos.R"); source("R/03-evaluacion.R")

nile <- leer_serie(Nile, fuente = "R datasets: Nile", unidad = "10^8 m^3")
part <- particionar(nile)                       # últimas h observaciones a validación
opt  <- optimizar(part$estimacion$y, "ses")     # alpha por rejilla, no a ojo
aj   <- ajustar_ses(part$estimacion$y, opt$optimo$alpha)

aj$yhat[1:5]            # pronósticos de un paso; NA en el calentamiento
aj$pronosticar(12)      # 12 pronósticos extramuestrales
aj$parametros$alpha     # parámetros y estados finales

e <- part$estimacion$y - aj$yhat
medidas(part$estimacion$y, aj$yhat, escala = escala_mase(part$estimacion$y))
validar_errores(e, p = 1)
```

## 4. Convenciones que fijan los números

- **ACF**: divisor único `T` en numerador y denominador (Definición 2.5), lo que
  la hace idéntica a `stats::acf()`; la diferencia máxima verificada es menor
  que `1e-12`.
- **Banda del correlograma**: `qnorm((1 + 0.95)/2)/sqrt(n)`, la línea `clim0` de
  `plot.acf`. No es un estadístico: es el intervalo asintótico de `r_h` bajo
  ruido blanco (Bartlett). En el correlograma de errores, `n` es el número de
  errores.
- **Inicializaciones**: media simple y SES, `Yhat_2 = Y_1`; media móvil, `k`
  períodos de calentamiento; doble media móvil, `2k-1`; Holt, `L_1 = Y_1` y
  `That_1 = 0`; tendencias, sin calentamiento.
- **Media simple recursiva**: `Yhat_{t+1} = (1/t) sum_{i<=t} Y_i`, no la media de
  toda la muestra.
- **Partición**: `h = min(12, floor(0.2 T))`, y al menos un ciclo en series
  estacionales.
- **MASE**: escalado con el MAE del ingenuo (o ingenuo estacional) de un paso
  dentro del tramo de estimación.
- **HAC**: núcleo de Bartlett con `floor(4 (T/100)^(2/9))` rezagos.
- **Tendencia exponencial**: `exp(a + theta t)` estima la mediana condicional;
  `corregir_sesgo = TRUE` multiplica por `exp(sigma2_ln/2)`.

## 5. Resumen de resultados

Los valores de esta tabla los produce `ejemplos.R` y quedan también en
`figs/tabla-resumen.csv`. Todos los parámetros salen de `optimizar()`, ninguno
se fijó a ojo. El referente es el pronóstico ingenuo, o el ingenuo estacional en
las series estacionales (`JohnsonJohnson` y `AirPassengers`).

| # | Serie | Método | Parámetros | h | MASE método | MASE referente |
|---|---|---|---|---|---|---|
| 1 | Nile | Media simple | ninguno | 12 | 0.843 | 0.835 |
| 2 | LakeHuron | Media móvil | k = 2 | 12 | 1.826 | 2.164 |
| 3 | Nile | SES | alpha = 0.24 | 12 | 0.806 | 0.835 |
| 4 | austres | Doble media móvil | k = 2 | 12 | 1.863 | 6.091 |
| 5 | austres | Tendencia lineal | MCO, 2 coef. | 12 | 3.974 | 6.091 |
| 6 | airmiles | Tendencia cuadrática | MCO, 3 coef. | 4 | 1.205 | 4.496 |
| 7 | JohnsonJohnson | Tendencia exponencial | MCO en logaritmos, con corrección de sesgo | 12 | 3.586 | 6.529 |
| 8 | WWWusage | Holt lineal | alpha = 0.95, beta = 0.95 | 12 | 1.890 | 7.625 |
| e | AirPassengers | Media simple (contraejemplo) | ninguno | 12 | 7.017 | 1.571 |

 En el ejemplo 1 la media simple **no supera** al ingenuo (0.843 contra 0.835): el
nivel del Nilo cambió tras 1899 y promediar toda la historia arrastra un nivel
antiguo. En el ejemplo 8 el óptimo de Holt cae en el borde de la rejilla
(alpha = beta = 0.95), lo que indica que el nivel y la pendiente de `WWWusage`
cambian casi tan rápido como la serie: el método se acerca a seguir el último
dato.

## 6. Declaración de uso de IA


- **Qué se pidió**: ayuda para estructurar el repositorio y para escribir una
  primera versión de las funciones `[...]`, así como la redacción de los
  encabezados de documentación.
- **Qué se recibió**: código en R para `[...]`, y el esqueleto de `informe.qmd`.
- **Qué se verificó por cuenta propia**: se contrastó la ACF propia contra
  `acf()` (diferencia máxima `[...]`), el estadístico de Ljung–Box contra
  `Box.test(type = "Ljung-Box")`, los coeficientes de las tendencias contra
  `lm()`, la equivalencia entre la forma de corrección de error y el promedio
  ponderado en SES, y las ecuaciones de Holt contra su forma de corrección de
  error. Además se reprodujeron a mano las cifras de las notas de la Clase 3
  sobre las series de las notas.
