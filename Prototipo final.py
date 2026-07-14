import streamlit as st
import pandas as pd
import numpy as np

# ==========================================
# CONFIGURACIÓN Y CONSTANTES
# ==========================================
# Columnas de resultados (del paciente y enzimas si existen)
COLUMNA_PACIENTE = 'AHG'  
COLUMNA_ENZIMA = 'ENZ'    

# Lista completa de antígenos del panel
ANTIGENOS_TODOS = [
    'D', 'C', 'E', 'c', 'e', 'Cw', 'M', 'N', 'S', 's', 'K', 'k', 
    'Kpa', 'Kpb', 'Jsa', 'Jsb', 'P1', 'Lea', 'Leb', 'Fya', 'Fyb', 
    'Jka', 'Jkb', 'Lua', 'Lub', 'Xga'
]

# Definición de parejas alélicas para evaluar dosis (Homocigosis vs Heterocigosis)
PAREJAS_CIGOTICAS = {
    'C': 'c', 'c': 'C', 'E': 'e', 'e': 'E',
    'M': 'N', 'N': 'M', 'S': 's', 's': 'S',
    'Fya': 'Fyb', 'Fyb': 'Fya', 'Jka': 'Jkb', 'Jkb': 'Jka',
    'K': 'k', 'k': 'K', 'Kpa': 'Kpb', 'Kpb': 'Kpa',
    'Jsa': 'Jsb', 'Jsb': 'Jsa', 'Lea': 'Leb', 'Leb': 'Lea',
    'Lua': 'Lub', 'Lub': 'Lua'
}

# Antígenos de alta frecuencia (se excluyen de la búsqueda de coincidencias completas estándar)
ALTA_FRECUENCIA = ['k', 'Kpb', 'Jsb', 'Lub']

# Antígenos de baja frecuencia (se descartan de las mezclas y búsquedas principales)
BAJA_FRECUENCIA = ['Cw', 'Kpa', 'Jsa', 'Lua']

# Efecto de las enzimas sobre los antígenos (D = Destruye, S = Sensibiliza/No afecta/Potencia)
EFECTO_ENZIMAS = {'Fya': 'D', 'Fyb': 'D', 'M': 'D', 'N': 'D', 'S': 'D', 's': 'D', 'Xga': 'D'}

# ==========================================
# INTERFAZ STREAMLIT
# ==========================================
st.title("Identificación de Anticuerpos Irregulares 🧪")
datos = None

archivo = st.file_uploader("Sube tu archivo CSV de panel", type=["csv"])

if archivo is not None:
    datos = pd.read_csv(archivo, delimiter=";")
    st.subheader("Vista previa de datos")
    st.dataframe(datos.head())

    # Limpieza
    columnas_criticas = [col for col in ANTIGENOS_TODOS+[COLUMNA_PACIENTE] if col in datos.columns]
    datos = datos.dropna(subset=columnas_criticas, how='all')
    datos = datos.rename(index={i:f"Célula {i+1}" for i in range(len(datos))})

    usar_enzimas = COLUMNA_ENZIMA in datos.columns and datos[COLUMNA_ENZIMA].notna().any()
    columnas_a_convertir = ANTIGENOS_TODOS+[COLUMNA_PACIENTE]
    if usar_enzimas:
        columnas_a_convertir.append(COLUMNA_ENZIMA)
    columnas_validas = [col for col in columnas_a_convertir if col in datos.columns]
    datos[columnas_validas] = datos[columnas_validas].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    resultados_paciente = datos[COLUMNA_PACIENTE]
    resultados_enzima = datos[COLUMNA_ENZIMA] if usar_enzimas else None
else:
    st.info("Por favor, sube un archivo CSV para iniciar el análisis. Ignora ese error feo de abajo, después lo arreglo.")

# ==========================================
#         FUNCIONES DE EVALUACIÓN
# ==========================================

def evaluar_dosis_mezcla(antigeno_evaluado, otro_sospechoso, df, resultados, col_intensidad):
    """
    Evalúa si la intensidad de reacción es mayor para homocigotas (dosis doble)
    que para heterocigotas, aislando el efecto del otro antígeno en la mezcla.
    Retorna la diferencia de medias (positiva si homocigota > heterocigota).
    """
    df_aislado = df[df[otro_sospechoso] == 0]
    resultados_aislados = resultados.loc[df_aislado.index]
    intensidades_aisladas = df_aislado.loc[resultados_aislados.index, col_intensidad]
    pareja = PAREJAS_CIGOTICAS.get(antigeno_evaluado)
    
    diff = 0.0
    if pareja and pareja in df_aislado.columns:
        # Homocigotas del antígeno evaluado: Antígeno=1 Y Pareja=0, y el resultado es positivo.
        mask_homo = (df_aislado[antigeno_evaluado] == 1) & (df_aislado[pareja] == 0)
        intensidades_homo = intensidades_aisladas[mask_homo]
        
        # Heterocigotas: Antígeno=1 Y Pareja=1, y el resultado es positivo.
        mask_hetero = (df_aislado[antigeno_evaluado] == 1) & (df_aislado[pareja] == 1)
        intensidades_hetero = intensidades_aisladas[mask_hetero]
        
        # Calcular diferencia si hay datos de ambos tipos
        if not intensidades_homo.empty and not intensidades_hetero.empty:
            diff = intensidades_homo.mean() - intensidades_hetero.mean()
            
    return diff


def validar_coherencia_dosis(ant1, ant2, df, resultados, col_intensidad):
    """
    Verifica si las reacciones del paciente con la mezcla tienen coherencia biológica.
    Si las intensidades del paciente son planas (todas iguales sin importar si es
    homocigoto o heterocigoto), la hipótesis de mezcla no es consistente (puede ser Anti-Alta Frecuencia).
    """
    diff1 = evaluar_dosis_mezcla(ant1, ant2, df, resultados, col_intensidad)
    diff2 = evaluar_dosis_mezcla(ant2, ant1, df, resultados, col_intensidad)
    
    # Si las diferencias de dosis para ambos antígenos son extremadamente bajas (cercanas a 0 o negativas)
    # y las intensidades del paciente son uniformemente reactivas, hay sospecha de anticuerpo único plano.
    casi_plano = abs(diff1) < 0.25 and abs(diff2) < 0.25
    return not casi_plano


def evaluar_alta_frecuencia(df, resultados, col_intensidad):
    """
    Busca si el perfil de resultados del paciente coincide con un anticuerpo
    dirigido contra un antígeno de alta frecuencia de forma única.
    """
    sospechosos_alta = []
    celulas_positivas = df[resultados > 0]
    celulas_negativas = df[resultados == 0]
    
    for ant in ALTA_FRECUENCIA:
        if ant in df.columns:
            # 1. Debe estar presente en todas (o casi todas) las células reactivas
            coincide_positivos = (celulas_positivas[ant] == 1).all()
            
            # 2. Si hay células negativas del paciente, el antígeno de alta frecuencia DEBE ser 0 en ellas (homocigosis negativa)
            coincide_negativos = True
            if not celulas_negativas.empty:
                coincide_negativos = (celulas_negativas[ant] == 0).all()
                
            if coincide_positivos and coincide_negativos:
                sospechosos_alta.append(ant)
                
    return sospechosos_alta


def imprimir_control_unico(antigeno, df, resultados):
    """
    Genera el texto simplificado para el control 3+3 de un anticuerpo único.
    """
    st.subheader("Control de confirmación 3+3 (Anticuerpo único):")
    pareja = PAREJAS_CIGOTICAS.get(antigeno)
    
    if pareja and pareja in df.columns:
        resultados_alineados = resultados.loc[df.index]
        n_homo_pos = len(df[(df[antigeno] == 1) & (df[pareja] == 0) & (resultados_alineados > 0)])
    else:
        resultados_alineados = resultados.loc[df.index]
        n_homo_pos = len(df[(df[antigeno] == 1) & (resultados_alineados > 0)])
        
    n_neg_no_reactivas = len(df[(df[antigeno] == 0) & (resultados_alineados == 0)])
    
    cumple = (n_homo_pos >= 3) and (n_neg_no_reactivas >= 3)
    estado = "Cumple" if cumple else "No cumple"
    
    st.write(f"[{estado}] Anti-{antigeno}: {n_homo_pos} células reactivas homocigotas (dosis doble) y {n_neg_no_reactivas} células negativas no reactivas.")


def imprimir_control_mezcla(antig_1, antig_2, df, resultados_paciente, col_ahg, col_enz, usar_enz):
    """
    Genera el control 3+3 adaptativo para mezclas. 
    Si uno se destruye y el otro se mantiene/potencia:
    - El resistente se confirma en la fase ENZ.
    - El destruido se confirma estrictamente en la fase basal (AHG) aislando al resistente.
    """
    st.subheader("Control de confirmación 3+3 (Mezcla de anticuerpos):")
    
    destruido_1 = EFECTO_ENZIMAS.get(antig_1) == 'D'
    destruido_2 = EFECTO_ENZIMAS.get(antig_2) == 'D'
    
    # Caso: Uno se destruye y el otro resiste (con columna ENZ disponible)
    if usar_enz and (destruido_1 != destruido_2):
        if destruido_1:
            destruido, resistente = antig_1, antig_2
        else:
            destruido, resistente = antig_2, antig_1
            
        pareja_res = PAREJAS_CIGOTICAS.get(resistente)
        pareja_dest = PAREJAS_CIGOTICAS.get(destruido)
        
        # --- 1. EVALUAR EL RESISTENTE EN ENZIMAS ---
        if pareja_res and pareja_res in df.columns:
            n_pos_res = len(df[(df[resistente] == 1) & (df[pareja_res] == 0) & (df[col_enz] > 0)])
        else:
            n_pos_res = len(df[(df[resistente] == 1) & (df[col_enz] > 0)])
            
        n_neg_res = len(df[(df[resistente] == 0) & (df[col_enz] == 0)])
        
        cumple_res = (n_pos_res >= 3) and (n_neg_res >= 3)
        est_res = "Cumple" if cumple_res else "No cumple"
        
        # --- 2. EVALUAR EL DESTRUIDO EN FASE BASAL (AHG) ---
        df_puro_dest = df[df[resistente] == 0]
        
        if pareja_dest and pareja_dest in df_puro_dest.columns:
            n_pos_dest = len(df_puro_dest[(df_puro_dest[destruido] == 1) & (df_puro_dest[pareja_dest] == 0) & (df_puro_dest[col_ahg] > 0)])
        else:
            n_pos_dest = len(df_puro_dest[(df_puro_dest[destruido] == 1) & (df_puro_dest[col_ahg] > 0)])
            
        n_neg_dest = len(df[(df[destruido] == 0) & (df[resistente] == 0) & (resultados_paciente == 0)])
        
        cumple_dest = (n_pos_dest >= 3) and (n_neg_dest >= 3)
        est_dest = "Cumple" if cumple_dest else "No cumple"
        
        st.write(f"[{est_res}] Anti-{resistente} (Evaluado en ENZ): {n_pos_res} células reactivas homocigotas y {n_neg_res} células negativas no reactivas.")
        st.write(f"[{est_dest}] Anti-{destruido} (Evaluado en AHG puro): {n_pos_dest} células reactivas homocigotas puras y {n_neg_dest} células negativas puras no reactivas.")

    else:
        # --- MODO TRADICIONAL (Si ambos resisten, se destruyen o no hay datos enzimáticos) ---
        pareja_1 = PAREJAS_CIGOTICAS.get(antig_1)
        df_puro_1 = df[df[antig_2] == 0]
        
        if pareja_1 and pareja_1 in df_puro_1.columns:
            n_homo_pos_1 = len(df_puro_1[(df_puro_1[antig_1] == 1) & (df_puro_1[pareja_1] == 0) & (df_puro_1[col_ahg] > 0)])
        else:
            n_homo_pos_1 = len(df_puro_1[(df_puro_1[antig_1] == 1) & (df_puro_1[col_ahg] > 0)])
            
        n_neg_no_reactivas_1 = len(df[(df[antig_1] == 0) & (df[antig_2] == 0) & (resultados_paciente == 0)])
        cumple_1 = (n_homo_pos_1 >= 3) and (n_neg_no_reactivas_1 >= 3)
        estado_1 = "Cumple" if cumple_1 else "No cumple"
        
        pareja_2 = PAREJAS_CIGOTICAS.get(antig_2)
        df_puro_2 = df[df[antig_1] == 0]
        
        if pareja_2 and pareja_2 in df_puro_2.columns:
           n_homo_pos_2 = len(df_puro_2[(df_puro_2[antig_2] == 1) & (df_puro_2[pareja_2] == 0) & (df_puro_2[col_ahg] > 0)])
        else:
            n_homo_pos_2 = len(df_puro_2[(df_puro_2[antig_2] == 1) & (df_puro_2[col_ahg] > 0)])
            
        n_neg_no_reactivas_2 = len(df[(df[antig_2] == 0) & (df[antig_1] == 0) & (resultados_paciente == 0)])
        cumple_2 = (n_homo_pos_2 >= 3) and (n_neg_no_reactivas_2 >= 3)
        estado_2 = "Cumple" if cumple_2 else "No cumple"
        
        st.write(f"[{estado_1}] Anti-{antig_1}: {n_homo_pos_1} células reactivas homocigotas puras y {n_neg_no_reactivas_1} células negativas puras no reactivas.")
        st.write(f"[{estado_2}] Anti-{antig_2}: {n_homo_pos_2} células reactivas homocigotas puras y {n_neg_no_reactivas_2} células negativas puras no reactivas.")



# ==========================================
#             EJECUCIÓN DEL ANÁLISIS
# ==========================================

# --- LIMPIEZA: Eliminar filas vacías y basura del CSV ---
columnas_criticas = [col for col in ANTIGENOS_TODOS + [COLUMNA_PACIENTE] if col in datos.columns]
datos = datos.dropna(subset=columnas_criticas, how='all')

# Renombrar las filas de las células activas
datos = datos.rename(index={i: f"Célula {i+1}" for i in range(len(datos))})

# 2. Pre-procesamiento de datos numéricos
usar_enzimas = COLUMNA_ENZIMA in datos.columns and datos[COLUMNA_ENZIMA].notna().any()

columnas_a_convertir = ANTIGENOS_TODOS + [COLUMNA_PACIENTE]
if usar_enzimas:
    columnas_a_convertir.append(COLUMNA_ENZIMA)

columnas_validas = [col for col in columnas_a_convertir if col in datos.columns]
datos[columnas_validas] = datos[columnas_validas].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

# Series de resultados limpios
resultados_paciente = datos[COLUMNA_PACIENTE]
resultados_enzima = datos[COLUMNA_ENZIMA] if usar_enzimas else None

# Variables para almacenar conclusiones finales
antig_confirmar_u = None
confirmar_mezcla = None

# ==========================================================
# NUEVA LÓGICA CLÍNICA: FILTRO DE DESCARTE POR NEGATIVOS
# ==========================================================

# 1. Identificar antígenos descartados (Presentes en células donde AHG == 0)
# Se excluyen de este descarte automático estricto a los de baja frecuencia
antigenos_descartados = set()
celulas_negativas = datos[resultados_paciente == 0]

for ant in ANTIGENOS_TODOS:
    if ant in datos.columns and ant not in BAJA_FRECUENCIA:
        # Si el antígeno está presente (1) en una célula donde la reacción real es 0, se descarta
        if (celulas_negativas[ant] == 1).any():
            antigenos_descartados.add(ant)

# 2. Identificar candidatos viables (Los que NO fueron descartados)
candidatos_no_descartados = [
    ant for ant in ANTIGENOS_TODOS 
    if ant in datos.columns and ant not in antigenos_descartados and ant not in ALTA_FRECUENCIA and ant not in BAJA_FRECUENCIA
]

# --- PASO 1: EVALUAR SI UN ÚNICO ANTICUERPO EXPLICA TODOS LOS POSITIVOS ---
coincidencias_completas = []
celulas_positivas = datos[resultados_paciente > 0]

for ant in candidatos_no_descartados:
    # Para ser considerado culpable único, debe estar presente en TODAS las células reactivas
    if (celulas_positivas[ant] == 1).all():
        coincidencias_completas.append(ant)

# Si hay exactamente un antígeno que sobrevivió al filtro y cubre los positivos, nos quedamos con él
if len(coincidencias_completas) == 1:
    antig_confirmar_u = coincidencias_completas[0]
    st.success(f"Resultado: Anti-{antig_confirmar_u}")

else:
    # --- PASO 2: SI NO HAY ÚNICO, INICIAR BÚSQUEDA DE MEZCLAS ---
    fallar_a_modo_basal = False
    
    if usar_enzimas:
        # --- MODO ENZIMAS (Mezclas asistidas por enzimas sobre los no descartados) ---
        celulas_negativizadas = datos[(resultados_paciente > 0) & (resultados_enzima == 0)]
        sospechosos_destruidos = []
        
        if not celulas_negativizadas.empty:
            for ant in candidatos_no_descartados:
                if EFECTO_ENZIMAS.get(ant) == 'D' and (celulas_negativizadas[ant] == 1).all():
                    sospechosos_destruidos.append(ant)
                        
        celulas_positivas_enz = datos[resultados_enzima > 0]
        sospechosos_resistentes = []
        if not celulas_positivas_enz.empty:
            for ant in candidatos_no_descartados:
                if EFECTO_ENZIMAS.get(ant) != 'D' and (celulas_positivas_enz[ant] == 1).all():
                    sospechosos_resistentes.append(ant)

        total_sospechosos = sospechosos_destruidos + sospechosos_resistentes
        
        if len(total_sospechosos) == 0:
            fallar_a_modo_basal = True
        elif len(total_sospechosos) == 1:
            antig_confirmar_u = total_sospechosos[0]
            st.success(f"Resultado: Anti-{antig_confirmar_u}")
        else:
            confirmar_mezcla = total_sospechosos
            # No imprimimos inmediatamente; validaremos la dosis en la sección posterior de impresión

    if not usar_enzimas or fallar_a_modo_basal:
        # --- MODO BASAL ESTÁNDAR (Puntuación probabilística sobre los candidatos restantes) ---
        evaluaciones_mezclas = []
        
        for i in range(len(candidatos_no_descartados)):
            for j in range(i+1, len(candidatos_no_descartados)):
                ant1 = candidatos_no_descartados[i]
                ant2 = candidatos_no_descartados[j]
                
                puntos_positivos_explicados = 0
                penalizaciones_falsos_positivos = 0
                
                for idx, fila in datos.iterrows():
                    reaccion_real = fila[COLUMNA_PACIENTE]
                    tiene_antigenos = (fila[ant1] == 1 or fila[ant2] == 1)
                    
                    if reaccion_real > 0 and tiene_antigenos:
                        puntos_positivos_explicados += 1
                    elif reaccion_real == 0 and tiene_antigenos:
                        penalizaciones_falsos_positivos += 1.5 
                
                score_cobertura = puntos_positivos_explicados - penalizaciones_falsos_positivos
                
                diff1 = evaluar_dosis_mezcla(ant1, ant2, datos, resultados_paciente, COLUMNA_PACIENTE)
                diff2 = evaluar_dosis_mezcla(ant2, ant1, datos, resultados_paciente, COLUMNA_PACIENTE)
                
                evaluaciones_mezclas.append({
                    'pareja': (ant1, ant2),
                    'score': score_cobertura,
                    'suma_diffs': diff1 + diff2
                })
        
        evaluaciones_mezclas = sorted(evaluaciones_mezclas, key=lambda x: (x['score'], x['suma_diffs']), reverse=True)
        
        if evaluaciones_mezclas and evaluaciones_mezclas[0]['score'] > 0:
            confirmar_mezcla = evaluaciones_mezclas[0]['pareja']
        else:
            # Resguardo: Si la puntuación de mezcla falla pero hay un antígeno sólido no descartado
            candidatos_validos_unicos = [ant for ant in candidatos_no_descartados if (celulas_positivas[ant] == 1).all()]
            if candidatos_validos_unicos:
                antig_confirmar_u = candidatos_validos_unicos[0]
            else:
                # Si fallan mezclas y no hay candidatos comunes, evaluamos si es uno de alta frecuencia
                sospechosos_alta = evaluar_alta_frecuencia(datos, resultados_paciente, COLUMNA_PACIENTE)
                if sospechosos_alta:
                    antig_confirmar_u = sospechosos_alta[0]
                    st.info("[SOPORTE ALTA FRECUENCIA] Detectado patrón compatible con anticuerpo de alta frecuencia.")
                else:
                    st.warning("Resultado: No se pudo determinar un anticuerpo o mezcla probable.")


# ==========================================================
# CONTROL DE COHERENCIA DE MEZCLA VS ALTA FRECUENCIA
# ==========================================================
if confirmar_mezcla and len(confirmar_mezcla) == 2:
    m1, m2 = confirmar_mezcla
    
    # Evaluar si la mezcla propuesta tiene coherencia en el efecto de dosis
    mezcla_es_coherente = validar_coherencia_dosis(m1, m2, datos, resultados_paciente, COLUMNA_PACIENTE)
    
    if not mezcla_es_coherente:
        # Si las intensidades son sospechosamente planas, descartamos la mezcla y evaluamos Alta Frecuencia
        st.warning(f"[FILTRO DE COHERENCIA] La mezcla propuesta (Anti-{m1} + Anti-{m2}) presenta un patrón de reacción plano (sin variación de dosis).")
        st.info("-> Re-evaluando sospecha de Anticuerpo contra Antígeno de Alta Frecuencia...")
        
        sospechosos_alta = evaluar_alta_frecuencia(datos, resultados_paciente, COLUMNA_PACIENTE)
        if sospechosos_alta:
            antig_confirmar_u = sospechosos_alta[0]
            confirmar_mezcla = None  # Cancelar la mezcla
            st.success(f"Resultado definitivo: Anti-{antig_confirmar_u}")
        else:
            # Si no se encontró un alta frecuencia claro, se mantiene la mezcla con una advertencia de patrón plano
            st.warning(f"Resultado (Mezcla más probable con advertencia): Anti-{m1} + Anti-{m2} (Reacciones inusualmente planas)")
    else:
        st.success(f"Resultado (Mezcla más probable): Anti-{m1} + Anti-{m2}")


# ==========================================
#         CONTROLES DE CONFIRMACIÓN 3+3
# ==========================================
if antig_confirmar_u:
    imprimir_control_unico(antig_confirmar_u, datos, resultados_paciente)
elif confirmar_mezcla and len(confirmar_mezcla) == 2:
    ant1, ant2 = confirmar_mezcla
    imprimir_control_mezcla(ant1, ant2, datos, resultados_paciente, COLUMNA_PACIENTE, COLUMNA_ENZIMA, usar_enzimas)
elif confirmar_mezcla and len(confirmar_mezcla) != 2:
    st.subheader("Control de confirmación 3+3 (Mezcla múltiple):")
    for susp in confirmar_mezcla:
        n_pos = len(datos[(datos[susp] == 1) & (resultados_paciente > 0)])
        n_neg_u = len(datos[(datos[confirmar_mezcla].sum(axis=1) == 0) & (resultados_paciente == 0)])
        
        cumple = (n_pos >= 3) and (n_neg_u >= 3)
        estado = "Cumple" if cumple else "No cumple"
        st.write(f"[{estado}] Anti-{susp}: {n_pos} células reactivas y {n_neg_u} células negativas puras no reactivas.")

