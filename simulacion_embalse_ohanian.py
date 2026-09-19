import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# MÓDULO 1: CÁLCULO DE VOLUMEN SEGÚN TESIS OHANIAN (CAPÍTULOS 3 Y 4)
# =============================================================================

class CalculadorVolumenPolyhedron:
    """
    Implementación directa del algoritmo de Proyección Central (Lien & Kajiya / Lin / Ohanian, Cap 3)
    y Slicing por Plano de Corte (Ohanian, Cap 4) para calcular el volumen de un recipiente 3D.
    """
    
    @staticmethod
    def determinante_3x3(pts):
        """Calcula el determinante Jacobian (Det(T)) de 3 puntos proyectados desde el origen (Capítulo 3)."""
        return np.linalg.det(pts)

    @classmethod
    def calcular_volumen_solido(cls, triangulos):
        """
        Calcula el volumen total de un poliedro cerrado sumando los tetraedros 
        firmados proyectados al origen: Volume = sum(|T| / 6) (Eq. 50 / Listing 1).
        """
        volumen_total = 0.0
        for tri in triangulos:
            # tri es una matriz (3, 3) con los 3 vértices del triángulo
            det_T = cls.determinante_3x3(tri)
            volumen_tetraedro = det_T / 6.0  # Firma automática por dirección del normal
            volumen_total += volumen_tetraedro
            
        return abs(volumen_total)

    @classmethod
    def recortar_y_tapar_malla(cls, triangulos, cota_agua):
        """
        Implementa el algoritmo de Slicing y Capping de Ohanian (Capítulo 4, Listing 3).
        Recorta la malla con el plano z = cota_agua y añade los 'cap facets' horizontales.
        """
        puntos_plano = np.array([0.0, 0.0, cota_agua])
        normal_plano = np.array([0.0, 0.0, -1.0])  # Orientado hacia abajo para conservar la parte inferior (agua)
        
        epsilon = 1e-8
        triangulos_recortados = []
        puntos_corte_tapa = []  # Para crear las facetas de tapa (capping)

        for tri in triangulos:
            # Clasificar los 3 vértices con el producto escalar q • n (Eq. 53, Fig. 17)
            dots = [np.dot(p - puntos_plano, normal_plano) for p in tri]
            
            # Caso 1: Todo el triángulo está bajo el agua (conservar entero)
            if all(d >= -epsilon for d in dots):
                triangulos_recortados.append(tri)
                
            # Caso 2: Todo el triángulo está sobre el agua (descartar)
            elif all(d <= epsilon for d in dots):
                continue
                
            # Caso 3: El triángulo cruza el plano del agua (Slicing)
            else:
                pts_bajo = []
                pts_interseccion = []
                
                n_pts = len(tri)
                for i in range(n_pts):
                    p1 = tri[i]
                    p2 = tri[(i + 1) % n_pts]
                    d1 = dots[i]
                    d2 = dots[(i + 1) % n_pts]
                    
                    if d1 >= -epsilon:
                        pts_bajo.append(p1)
                        
                    # Si el segmento cruza el plano, calcular intersección vectorialmente (Eq. 53)
                    if (d1 > epsilon and d2 < -epsilon) or (d1 < -epsilon and d2 > epsilon):
                        r = p2 - p1
                        q = p1 - puntos_plano
                        # Factor de escala por proyecciones (Eq. 53)
                        ratio = abs(np.dot(q, normal_plano) / np.dot(r, normal_plano))
                        p_int = p1 + ratio * r
                        
                        pts_bajo.append(p_int)
                        pts_interseccion.append(p_int)
                        puntos_corte_tapa.append(p_int)

                # Triangular la polígonos resultantes bajo el agua
                if len(pts_bajo) >= 3:
                    for k in range(1, len(pts_bajo) - 1):
                        triangulos_recortados.append([pts_bajo[0], pts_bajo[k], pts_bajo[k+1]])

        # Algoritmo de Capping Facets (Tesis Cap. 4, Fig. 19, Table 3):
        # Tapa el poliedro abierto en el plano del agua usando abanicos de triángulos
        if len(puntos_corte_tapa) >= 3:
            # Centroide del espejo de agua para triangulación estable
            centroide_tapa = np.mean(puntos_corte_tapa, axis=0)
            centroide_tapa[2] = cota_agua  # Asegurar cota exactitud
            
            # Ordenar puntos angularmente alrededor del centroide
            puntos_2d = [p[:2] - centroide_tapa[:2] for p in puntos_corte_tapa]
            angulos = [np.arctan2(p[1], p[0]) for p in puntos_2d]
            indices_ordenados = np.argsort(angulos)
            pts_tapa_ord = [puntos_corte_tapa[idx] for idx in indices_ordenados]
            
            # Crear facetas de tapa en sentido correcto (outward normal)
            for i in range(len(pts_tapa_ord)):
                p_actual = pts_tapa_ord[i]
                p_siguiente = pts_tapa_ord[(i + 1) % len(pts_tapa_ord)]
                # Cap facet: {Centroide, Siguiente, Actual}
                triangulos_recortados.append([centroide_tapa, p_siguiente, p_actual])

        return triangulos_recortados


# =============================================================================
# MÓDULO 2: MODELADO DEL EMBALSE Y SIMULACIÓN HIDROLÓGICA
# =============================================================================

class EmbalseOhanian3D:
    def __init__(self, triangulos_malla_3d: list, cota_corona_vertedero: float):
        """
        :param triangulos_malla_3d: Lista de triángulos 3D (cada uno [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]])
        :param cota_corona_vertedero: Cota donde el embalse está al 100% de capacidad
        """
        self.malla_original = [np.array(t) for t in triangulos_malla_3d]
        self.cota_max = cota_corona_vertedero
        
        # Encontrar el fondo del embalse (z_min)
        todos_puntos = np.vstack([t for t in self.malla_original])
        self.z_min = np.min(todos_puntos[:, 2])
        
        # Capacidad total al 100%
        malla_full = CalculadorVolumenPolyhedron.recortar_y_tapar_malla(self.malla_original, self.cota_max)
        self.capacidad_100_m3 = CalculadorVolumenPolyhedron.calcular_volumen_solido(malla_full)

    def obtener_volumen_en_cota(self, cota_agua: float) -> float:
        """Devuelve el volumen almacenado a una cota específica."""
        if cota_agua <= self.z_min:
            return 0.0
        malla_recortada = CalculadorVolumenPolyhedron.recortar_y_tapar_malla(self.malla_original, cota_agua)
        return CalculadorVolumenPolyhedron.calcular_volumen_solido(malla_recortada)

    def simular_llenado_por_lluvia(self, 
                                   area_cuenca_km2: float, 
                                   c_escorrentia: float, 
                                   precipitacion_mm_h: float, 
                                   paso_minutos: float = 10.0):
        """
        Simula el llenado por evento de lluvia hasta sobrepasar el 100% de capacidad.
        """
        # Caudal aportante (Método Racional: Q = C * I * A)
        area_m2 = area_cuenca_km2 * 1_000_000.0
        intensidad_m_s = (precipitacion_mm_h / 1000.0) / 3600.0
        caudal_m3s = c_escorrentia * intensidad_m_s * area_m2
        
        dt_seg = paso_minutos * 60.0
        volumen_actual = 0.0
        tiempo_acumulado_min = 0.0
        
        # Discretizar curva Elevación-Volumen
        cotas_curva = np.linspace(self.z_min, self.cota_max + 1.0, 30)
        vols_curva = [self.obtener_volumen_en_cota(c) for c in cotas_curva]
        
        registros = []
        
        print("=" * 65)
        print("   SIMULADOR DE EMBALSE BASADO EN TESIS OHANIAN (VIRGINIA TECH)")
        print("=" * 65)
        print(f"Capacidad Máxima (100%):   {self.capacidad_100_m3:,.2f} m³")
        print(f"Cota Vertedero (100%):     {self.cota_max:.2f} m.s.n.m.")
        print(f"Precipitación Aplicada:     {precipitacion_mm_h:.1f} mm/h")
        print(f"Caudal Aportante Entrante:  {caudal_m3s:.2f} m³/s ({caudal_m3s*3600:,.1f} m³/h)")
        print("-" * 65)
        
        tiempo_100 = None
        
        while True:
            porcentaje = (volumen_actual / self.capacidad_100_m3) * 100.0
            cota_actual = np.interp(volumen_actual, vols_curva, cotas_curva)
            
            registros.append({
                'tiempo_min': tiempo_acumulado_min,
                'tiempo_horas': tiempo_acumulado_min / 60.0,
                'volumen_m3': volumen_actual,
                'porcentaje_%': porcentaje,
                'cota_m': cota_actual
            })
            
            if volumen_actual >= self.capacidad_100_m3 and tiempo_100 is None:
                tiempo_100 = tiempo_acumulado_min
                print(f"\n ALERTA: Embalse al 100% de capacidad superado a las {tiempo_acumulado_min/60.0:.2f} horas "
                      f"({tiempo_acumulado_min:.0f} minutos).")
                break
                
            volumen_actual += caudal_m3s * dt_seg
            tiempo_acumulado_min += paso_minutos

        df_res = pd.DataFrame(registros)
        return df_res, tiempo_100


# =============================================================================
# PRUEBA Y DEMOSTRACIÓN (CONSTRUCCIÓN DE MALLA 3D SINTÉTICA DE VALLE)
# =============================================================================
if __name__ == "__main__":

    # 1. GENERACIÓN DE UNA MALLA POLIEDRAL 3D DEL EMBALSE (Forma cónica/valle)
    def generar_malla_embalse_3d(r_max=100.0, h_max=30.0, cota_base=100.0, n_div=24):
        triangulos = []
        angles = np.linspace(0, 2*np.pi, n_div, endpoint=False)
        
        # Puntos base (fondo) y superiores (borde)
        for i in range(n_div):
            a1 = angles[i]
            a2 = angles[(i + 1) % n_div]
            
            # Vértices del fondo (cota 100)
            p_fondo_1 = [0.1 * r_max * np.cos(a1), 0.1 * r_max * np.sin(a1), cota_base]
            p_fondo_2 = [0.1 * r_max * np.cos(a2), 0.1 * r_max * np.sin(a2), cota_base]
            
            # Vértices de la pared alta (cota 130)
            p_alto_1 = [r_max * np.cos(a1), r_max * np.sin(a1), cota_base + h_max]
            p_alto_2 = [r_max * np.cos(a2), r_max * np.sin(a2), cota_base + h_max]
            
            # Triángulo de fondo (suelo)
            triangulos.append([[0,0,cota_base], p_fondo_2, p_fondo_1])
            # Paredes laterales
            triangulos.append([p_fondo_1, p_fondo_2, p_alto_2])
            triangulos.append([p_fondo_1, p_alto_2, p_alto_1])
            
        return triangulos

    # Crear malla 3D del vaso
    malla_3d = generar_malla_embalse_3d()
    COTA_DESBORDE = 125.0  # Cota máxima (100% de capacidad)

    # 2. INSTANCIAR EMBALSE CON LA METODOLOGÍA OHANIAN
    embalse = EmbalseOhanian3D(malla_3d, cota_corona_vertedero=COTA_DESBORDE)

    # 3. EJECUTAR SIMULACIÓN CON PRECIPITACIÓN
    df_sim, tiempo_desborde = embalse.simular_llenado_por_lluvia(
        area_cuenca_km2=12.0,      # Cuenca de 12 km²
        c_escorrentia=0.50,       # Coeficiente de escorrentía
        precipitacion_mm_h=20.0,  # Lluvia de 20 mm/h
        paso_minutos=15.0         # Evaluación cada 15 min
    )

    # 4. GRAFICAR RESULTADOS DE LLENADO Y CAPACIDAD
    fig, ax1 = plt.subplots(figsize=(10, 5))

    color = 'tab:blue'
    ax1.set_xlabel('Tiempo de Lluvia (Horas)', fontweight='bold')
    ax1.set_ylabel('Porcentaje de Capacidad (%)', color=color, fontweight='bold')
    ax1.plot(df_sim['tiempo_horas'], df_sim['porcentaje_%'], color=color, linewidth=2.5, label='Llenado (%)')
    ax1.axhline(100, color='red', linestyle='--', linewidth=1.8, label='Capacidad Máxima (100%)')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax1.twinx()
    color = 'tab:green'
    ax2.set_ylabel('Cota del Agua (m.s.n.m.)', color=color, fontweight='bold')
    ax2.plot(df_sim['tiempo_horas'], df_sim['cota_m'], color=color, linestyle=':', linewidth=2, label='Cota del Agua')
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title('Simulación de Llenado de Embalse (Basado en Tesis Ohanian, VTech)', fontweight='bold')
    fig.tight_layout()
    plt.show()
