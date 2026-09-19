# Simulación de Embalse — Método de Malla 3D (Tesis Ohanian)

Modelo en Python para calcular el volumen almacenado en un embalse a partir
de una malla poliedral 3D (triángulos), y simular el tiempo que tardaría en
desbordarse ante un evento de lluvia, usando el **Método Racional**.

Basado en los algoritmos de la tesis de Ohanian (Virginia Tech) sobre cálculo
de volumen de sólidos mediante proyección central y *slicing* con capping de
malla (Capítulos 3 y 4).

## ¿Qué hace?

1. **`CalculadorVolumenPolyhedron`**:
   - Calcula el volumen de un sólido cerrado sumando tetraedros firmados
     proyectados al origen (método del determinante 3x3).
   - Recorta (*slicing*) la malla 3D con el plano de la superficie del agua
     y "tapa" el corte (*capping*) generando triángulos nuevos para cerrar
     el sólido resultante.

2. **`EmbalseOhanian3D`**: recibe la malla 3D del vaso del embalse y calcula,
   para cualquier cota de agua, el volumen almacenado usando el recorte +
   capping anterior.

3. **Simulación de lluvia** (`simular_llenado_por_lluvia`): aplica el Método
   Racional (Q = C·I·A) para estimar el caudal de entrada, integra el
   volumen en el tiempo, interpola la cota correspondiente y detecta el
   momento en que el embalse alcanza el 100% de su capacidad.

4. **Visualización**: gráfico con la evolución del porcentaje de llenado y
   la cota del agua a lo largo del tiempo.

## Requisitos

- Python 3.8+
- numpy
- pandas
- matplotlib

Instalación:

```bash
pip install -r requirements.txt
```

## Uso

```bash
python simulacion_embalse_ohanian.py
```

El script genera automáticamente una malla 3D sintética con forma de valle
(cónica), instancia el embalse, corre una simulación con una cuenca de
12 km², coeficiente de escorrentía de 0.50 y lluvia de 20 mm/h, y muestra
el resultado en consola junto con una gráfica.

## Personalización

En el bloque `if __name__ == "__main__":` puedes ajustar:

- `generar_malla_embalse_3d()`: parámetros de la malla sintética (`r_max`,
  `h_max`, `cota_base`, `n_div`), o reemplazarla por una malla real
  (por ejemplo, exportada desde un modelo de terreno / TIN).
- `COTA_DESBORDE`: cota del vertedero (100% de capacidad).
- Parámetros de la tormenta: `area_cuenca_km2`, `c_escorrentia`,
  `precipitacion_mm_h`, `paso_minutos`.

## Licencia

MIT
