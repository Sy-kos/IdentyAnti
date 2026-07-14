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
    pareja = PAREJAS_CIGOTICAS.get(antigeno)
    resultados_alineados = resultados.loc[df.index]
    if pareja and pareja in df.columns:
        n_homo_pos = len(df[(df[antigeno]==1)&(df[pareja]==0)&(resultados_alineados>0)])
    else:
        n_homo_pos = len(df[(df[antigeno]==1)&(resultados_alineados>0)])
    n_neg_no_reactivas = len(df[(df[antigeno]==0)&(resultados_alineados==0)])
    cumple = (n_homo_pos>=3) and (n_neg_no_reactivas>=3)
    estado = "Cumple" if cumple else "No cumple"
    return [f"[{estado}] Anti-{antigeno}: {n_homo_pos} células reactivas homocigotas y {n_neg_no_reactivas} negativas no reactivas."]


def imprimir_control_mezcla(antig_1, antig_2, df, resultados_paciente, col_ahg, col_enz, usar_enz):
    """
    Genera el control 3+3 adaptativo para mezclas.
    Retorna una lista de strings con los resultados en lugar de imprimir directamente.
    """
    salida = []
    salida.append("Control de confirmación 3+3 (Mezcla de anticuerpos):")

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
        est_res = "Cumple" if (n_pos_res >= 3 and n_neg_res >= 3) else "No cumple"

        # --- 2. EVALUAR EL DESTRUIDO EN FASE BASAL (AHG) ---
        df_puro_dest = df[df[resistente] == 0]
        if pareja_dest and pareja_dest in df_puro_dest.columns:
            n_pos_dest = len(df_puro_dest[(df_puro_dest[destruido] == 1) & (df_puro_dest[pareja_dest] == 0) & (df_puro_dest[col_ahg] > 0)])
        else:
            n_pos_dest = len(df_puro_dest[(df_puro_dest[destruido] == 1) & (df_puro_dest[col_ahg] > 0)])

        n_neg_dest = len(df[(df[destruido] == 0) & (df[resistente] == 0) & (resultados_paciente == 0)])
        est_dest = "Cumple" if (n_pos_dest >= 3 and n_neg_dest >= 3) else "No cumple"

        salida.append(f"[{est_res}] Anti-{resistente} (Evaluado en ENZ): {n_pos_res} células reactivas homocigotas y {n_neg_res} células negativas no reactivas.")
        salida.append(f"[{est_dest}] Anti-{destruido} (Evaluado en AHG puro): {n_pos_dest} células reactivas homocigotas puras y {n_neg_dest} células negativas puras no reactivas.")

    else:
        # --- MODO TRADICIONAL ---
        pareja_1 = PAREJAS_CIGOTICAS.get(antig_1)
        df_puro_1 = df[df[antig_2] == 0]
        if pareja_1 and pareja_1 in df_puro_1.columns:
            n_homo_pos_1 = len(df_puro_1[(df_puro_1[antig_1] == 1) & (df_puro_1[pareja_1] == 0) & (df_puro_1[col_ahg] > 0)])
        else:
            n_homo_pos_1 = len(df_puro_1[(df_puro_1[antig_1] == 1) & (df_puro_1[col_ahg] > 0)])

        n_neg_no_reactivas_1 = len(df[(df[antig_1] == 0) & (df[antig_2] == 0) & (resultados_paciente == 0)])
        estado_1 = "Cumple" if (n_homo_pos_1 >= 3 and n_neg_no_reactivas_1 >= 3) else "No cumple"

        pareja_2 = PAREJAS_CIGOTICAS.get(antig_2)
        df_puro_2 = df[df[antig_1] == 0]
        if pareja_2 and pareja_2 in df_puro_2.columns:
            n_homo_pos_2 = len(df_puro_2[(df_puro_2[antig_2] == 1) & (df_puro_2[pareja_2] == 0) & (df_puro_2[col_ahg] > 0)])
        else:
            n_homo_pos_2 = len(df_puro_2[(df_puro_2[antig_2] == 1) & (df_puro_2[col_ahg] > 0)])

        n_neg_no_reactivas_2 = len(df[(df[antig_2] == 0) & (df[antig_1] == 0) & (resultados_paciente == 0)])
        estado_2 = "Cumple" if (n_homo_pos_2 >= 3 and n_neg_no_reactivas_2 >= 3) else "No cumple"

        salida.append(f"[{estado_1}] Anti-{antig_1}: {n_homo_pos_1} células reactivas homocigotas puras y {n_neg_no_reactivas_1} células negativas puras no reactivas.")
        salida.append(f"[{estado_2}] Anti-{antig_2}: {n_homo_pos_2} células reactivas homocigotas puras y {n_neg_no_reactivas_2} células negativas puras no reactivas.")

    return salida

# ==========================================
# INTERFAZ STREAMLIT
# ==========================================
st.title("Identificación de Anticuerpos Irregulares 🧪")

antig_confirmar_u = None
confirmar_mezcla = None
conclusion = None
controles = []


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

    # ============================
    # EJECUCIÓN DE LA LÓGICA CLÍNICA
    # ============================
    antig_confirmar_u = None
    confirmar_mezcla = None
    conclusion = None
    controles = []

    # Filtro de descarte
    antigenos_descartados = set()
    celulas_negativas = datos[resultados_paciente == 0]
    for ant in ANTIGENOS_TODOS:
        if ant in datos.columns and ant not in BAJA_FRECUENCIA:
            if (celulas_negativas[ant] == 1).any():
                antigenos_descartados.add(ant)

    candidatos_no_descartados = [
        ant for ant in ANTIGENOS_TODOS
        if ant in datos.columns and ant not in antigenos_descartados
        and ant not in ALTA_FRECUENCIA and ant not in BAJA_FRECUENCIA
    ]

    # Paso 1: único anticuerpo
    coincidencias_completas = []
    celulas_positivas = datos[resultados_paciente > 0]
    for ant in candidatos_no_descartados:
        if (celulas_positivas[ant] == 1).all():
            coincidencias_completas.append(ant)

    if len(coincidencias_completas) == 1:
        antig_confirmar_u = coincidencias_completas[0]
        conclusion = f"Resultado: Anti-{antig_confirmar_u}"
    else:
        # Paso 2: mezclas con enzimas
        sospechosos_destruidos = []
        sospechosos_resistentes = []
        if usar_enzimas:
            # Células que eran positivas en AHG pero se negativizaron en ENZ
            celulas_negativizadas = datos[(resultados_paciente > 0) & (resultados_enzima == 0)]
            if not celulas_negativizadas.empty:
                for ant in candidatos_no_descartados:
                    if EFECTO_ENZIMAS.get(ant) == 'D' and (celulas_negativizadas[ant] == 1).all():
                        sospechosos_destruidos.append(ant)

            # Células que siguen positivas en ENZ
            celulas_positivas_enz = datos[resultados_enzima > 0]
            if not celulas_positivas_enz.empty:
                for ant in candidatos_no_descartados:
                    if EFECTO_ENZIMAS.get(ant) != 'D' and (celulas_positivas_enz[ant] == 1).all():
                        sospechosos_resistentes.append(ant)

        # Unimos los sospechosos
        total_sospechosos = sospechosos_destruidos + sospechosos_resistentes

        if len(total_sospechosos) == 1:
            antig_confirmar_u = total_sospechosos[0]
            conclusion = f"Resultado: Anti-{antig_confirmar_u}"
        elif len(total_sospechosos) == 2:
            confirmar_mezcla = total_sospechosos

    # Si no hay sospechosos claros, pasamos a evaluar alta frecuencia
    sospechosos_alta = evaluar_alta_frecuencia(datos, resultados_paciente, COLUMNA_PACIENTE)
    if sospechosos_alta:
        antig_confirmar_u = sospechosos_alta[0]
        conclusion = f"[Advertencia antígeno de alta frecuencia] Anti-{antig_confirmar_u}"
    else:
        # --- MODO BASAL ESTÁNDAR ---
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
            # Resguardo: intentar un candidato único
            candidatos_validos_unicos = [
                ant for ant in candidatos_no_descartados
                if (celulas_positivas[ant] == 1).all()
            ]
            if candidatos_validos_unicos:
                antig_confirmar_u = candidatos_validos_unicos[0]
                conclusion = f"Resultado: Anti-{antig_confirmar_u}"
            else:
                sospechosos_alta = evaluar_alta_frecuencia(datos, resultados_paciente, COLUMNA_PACIENTE)
                if sospechosos_alta:
                    antig_confirmar_u = sospechosos_alta[0]
                    conclusion = f"[SOPORTE ALTA FRECUENCIA] Anti-{antig_confirmar_u}"
                else:
                    conclusion = "Resultado: No se pudo determinar un anticuerpo o mezcla probable."# ============================
 # ============================
 # VALIDACIÓN DE MEZCLA
 # ============================
 if confirmar_mezcla is not None and len(confirmar_mezcla) == 2 and antig_confirmar_u is None:
    m1, m2 = confirmar_mezcla
    mezcla_es_coherente = validar_coherencia_dosis(m1, m2, datos, resultados_paciente, COLUMNA_PACIENTE)
    if mezcla_es_coherente:
        conclusion = f"Resultado (Mezcla más probable): Anti-{m1} + Anti-{m2}"
    else:
        sospechosos_alta = evaluar_alta_frecuencia(datos, resultados_paciente, COLUMNA_PACIENTE)
        if sospechosos_alta:
            antig_confirmar_u = sospechosos_alta[0]
            conclusion = f"Resultado definitivo: Anti-{antig_confirmar_u}"
        else:
            conclusion = f"Resultado (Mezcla con advertencia): Anti-{m1} + Anti-{m2} (patrón plano)"

 # ============================
 # CONTROLES DE CONFIRMACIÓN 3+3
 # ============================
 if antig_confirmar_u:
    resultado_unico = imprimir_control_unico(antig_confirmar_u, datos, resultados_paciente)
    controles.extend(resultado_unico)
    # Ajustar la conclusión según el control
    if "[Cumple]" in resultado_unico[0]:
        conclusion = f"Resultado definitivo: Anti-{antig_confirmar_u}"
    else:
        conclusion = f"Resultado tentativo: Anti-{antig_confirmar_u} (no cumple 3+3)"

 elif confirmar_mezcla and len(confirmar_mezcla) == 2:
    ant1, ant2 = confirmar_mezcla
    resultado_mezcla = imprimir_control_mezcla(ant1, ant2, datos, resultados_paciente, COLUMNA_PACIENTE, COLUMNA_ENZIMA, usar_enzimas)
    controles.extend(resultado_mezcla)
    # Ajustar conclusión según controles
    if any("[Cumple]" in r for r in resultado_mezcla):
        conclusion = f"Resultado (Mezcla más probable): Anti-{ant1} + Anti-{ant2}"
    else:
        conclusion = f"Resultado tentativo: Anti-{ant1} + Anti-{ant2} (no cumple 3+3)"

 elif confirmar_mezcla and len(confirmar_mezcla) > 2:
    for susp in confirmar_mezcla:
        n_pos = len(datos[(datos[susp] == 1) & (resultados_paciente > 0)])
        n_neg_u = len(datos[(datos[confirmar_mezcla].sum(axis=1) == 0) & (resultados_paciente == 0)])
        cumple = (n_pos >= 3) and (n_neg_u >= 3)
        estado = "Cumple" if cumple else "No cumple"
        controles.append(f"[{estado}] Anti-{susp}: {n_pos} células reactivas y {n_neg_u} negativas puras no reactivas")
    # Ajustar conclusión según si alguno cumple
    if any("Cumple" in c for c in controles):
        conclusion = f"Resultado (Mezcla probable): {', '.join(confirmar_mezcla)}"
    else:
        conclusion = f"Resultado tentativo: {', '.join(confirmar_mezcla)} (no cumple 3+3)"

# ============================
# SALIDA EN STREAMLIT
# ============================
st.subheader("Conclusión")
if conclusion:
    st.success(conclusion)
else:
    st.warning("No se pudo determinar un resultado.")

if controles:
    st.subheader("Controles de confirmación 3+3")
    for c in controles:
        st.write(c)
