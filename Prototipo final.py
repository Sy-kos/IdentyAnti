import streamlit as st
import pandas as pd
import numpy as np

# ==========================================
# CONFIGURACIÓN Y CONSTANTES
# ==========================================
COLUMNA_PACIENTE = 'AHG'
COLUMNA_ENZIMA = 'ENZ'

ANTIGENOS_TODOS = [
    'D','C','E','c','e','Cw','M','N','S','s','K','k',
    'Kpa','Kpb','Jsa','Jsb','P1','Lea','Leb','Fya','Fyb',
    'Jka','Jkb','Lua','Lub','Xga'
]

PAREJAS_CIGOTICAS = {
    'C':'c','c':'C','E':'e','e':'E','M':'N','N':'M','S':'s','s':'S',
    'Fya':'Fyb','Fyb':'Fya','Jka':'Jkb','Jkb':'Jka','K':'k','k':'K',
    'Kpa':'Kpb','Kpb':'Kpa','Jsa':'Jsb','Jsb':'Jsa','Lea':'Leb','Leb':'Lea',
    'Lua':'Lub','Lub':'Lua'
}

ALTA_FRECUENCIA = ['k','Kpb','Jsb','Lub']
BAJA_FRECUENCIA = ['Cw','Kpa','Jsa','Lua']
EFECTO_ENZIMAS = {'Fya':'D','Fyb':'D','M':'D','N':'D','S':'D','s':'D','Xga':'D'}

# ==========================================
# FUNCIONES DE EVALUACIÓN
# ==========================================
def evaluar_dosis_mezcla(antigeno_evaluado, otro_sospechoso, df, resultados, col_intensidad):
    df_aislado = df[df[otro_sospechoso] == 0]
    resultados_aislados = resultados.loc[df_aislado.index]
    intensidades_aisladas = df_aislado.loc[resultados_aislados.index, col_intensidad]
    pareja = PAREJAS_CIGOTICAS.get(antigeno_evaluado)
    diff = 0.0
    if pareja and pareja in df_aislado.columns:
        mask_homo = (df_aislado[antigeno_evaluado]==1)&(df_aislado[pareja]==0)
        intensidades_homo = intensidades_aisladas[mask_homo]
        mask_hetero = (df_aislado[antigeno_evaluado]==1)&(df_aislado[pareja]==1)
        intensidades_hetero = intensidades_aisladas[mask_hetero]
        if not intensidades_homo.empty and not intensidades_hetero.empty:
            diff = intensidades_homo.mean() - intensidades_hetero.mean()
    return diff

def validar_coherencia_dosis(ant1, ant2, df, resultados, col_intensidad):
    diff1 = evaluar_dosis_mezcla(ant1, ant2, df, resultados, col_intensidad)
    diff2 = evaluar_dosis_mezcla(ant2, ant1, df, resultados, col_intensidad)
    casi_plano = abs(diff1)<0.25 and abs(diff2)<0.25
    return not casi_plano

def evaluar_alta_frecuencia(df, resultados, col_intensidad):
    sospechosos_alta = []
    celulas_positivas = df[resultados>0]
    celulas_negativas = df[resultados==0]
    for ant in ALTA_FRECUENCIA:
        if ant in df.columns:
            coincide_positivos = (celulas_positivas[ant]==1).all()
            coincide_negativos = True
            if not celulas_negativas.empty:
                coincide_negativos = (celulas_negativas[ant]==0).all()
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
    return f"[{estado}] Anti-{antigeno}: {n_homo_pos} células reactivas homocigotas y {n_neg_no_reactivas} negativas no reactivas."

def imprimir_control_mezcla(antig_1, antig_2, df, resultados_paciente, col_ahg, col_enz, usar_enz):
    salida = []
    destruido_1 = EFECTO_ENZIMAS.get(antig_1) == 'D'
    destruido_2 = EFECTO_ENZIMAS.get(antig_2) == 'D'
    if usar_enz and (destruido_1 != destruido_2):
        if destruido_1:
            destruido, resistente = antig_1, antig_2
        else:
            destruido, resistente = antig_2, antig_1
        salida.append(f"[Mezcla con enzimas] Anti-{resistente} confirmado en ENZ, Anti-{destruido} confirmado en AHG")
    else:
        salida.append(f"[Mezcla tradicional] Anti-{antig_1} + Anti-{antig_2}")
    return salida

# ==========================================
# INTERFAZ STREAMLIT
# ==========================================
st.title("Identificación de Anticuerpos Irregulares 🧪")

opcion = st.radio("¿Qué quieres subir?", ["CSV", "Imagen"])

if opcion == "CSV":
    archivo = st.file_uploader("Sube tu archivo CSV de panel", type=["csv"])
    if archivo is not None:
        datos = pd.read_csv(archivo, delimiter=";")
        st.subheader("Vista previa de datos")
        st.dataframe(datos.head())

elif opcion == "Imagen":
    st.warning("Todavía no se con imagen, se hace lo que se puede.")
    imagen = st.file_uploader("Sube una imagen del panel", type=["png","jpg","jpeg"])
    if imagen is not None:
        # Aquí podrías dejarlo vacío o simular un DataFrame
        # datos = convertir_imagen_a_dataframe(imagen)
        st.info("Todavía mi cerebro no da para poner con imagen, sabrá dios cómo logré lo que hay, amén. Calmate por deoooos, siempre puedes modificar el csv manual, no seas floj@")


if 'datos' in locals():
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
                    conclusion = "Resultado: No se pudo determinar un anticuerpo o mezcla probable."
 
    # ============================
    # VALIDACIÓN DE MEZCLA
    # ============================
    if confirmar_mezcla and len(confirmar_mezcla) == 2:
        m1, m2 = confirmar_mezcla
        mezcla_es_coherente = validar_coherencia_dosis(m1, m2, datos, resultados_paciente, COLUMNA_PACIENTE)
        if mezcla_es_coherente:
            conclusion = f"Resultado (Mezcla más probable): Anti-{m1} + Anti-{m2}"
            controles.extend(
                imprimir_control_mezcla(
                    m1, m2,
                    datos, resultados_paciente,
                    COLUMNA_PACIENTE, COLUMNA_ENZIMA,
                    usar_enzimas
                )
            )
        else:
            sospechosos_alta = evaluar_alta_frecuencia(datos, resultados_paciente, COLUMNA_PACIENTE)
            if sospechosos_alta:
                antig_confirmar_u = sospechosos_alta[0]
                conclusion = f"Resultado definitivo: Anti-{antig_confirmar_u}"
                controles.append(imprimir_control_unico(antig_confirmar_u, datos, resultados_paciente))
            else:
                conclusion = f"Resultado (Mezcla con advertencia): Anti-{m1} + Anti-{m2} (patrón plano)"

    # ============================
    # CONTROLES DE CONFIRMACIÓN 3+3
    # ============================
    if antig_confirmar_u:
        controles.append(imprimir_control_unico(antig_confirmar_u, datos, resultados_paciente))
    elif confirmar_mezcla and len(confirmar_mezcla) == 2:
        ant1, ant2 = confirmar_mezcla
        controles.extend(imprimir_control_mezcla(ant1, ant2, datos, resultados_paciente, COLUMNA_PACIENTE, COLUMNA_ENZIMA, usar_enzimas))
    elif confirmar_mezcla and len(confirmar_mezcla) > 2:
        for susp in confirmar_mezcla:
            n_pos = len(datos[(datos[susp] == 1) & (resultados_paciente > 0)])
            n_neg_u = len(datos[(datos[confirmar_mezcla].sum(axis=1) == 0) & (resultados_paciente == 0)])
            cumple = (n_pos >= 3) and (n_neg_u >= 3)
            estado = "Cumple" if cumple else "No cumple"
            controles.append(f"[{estado}] Anti-{susp}: {n_pos} células reactivas y {n_neg_u} negativas puras no reactivas")

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

