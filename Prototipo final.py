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
