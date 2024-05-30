import pandas as pd
import numpy as np
import math
from sympy import symbols, diff
from sympy.solvers import solve



"""
2.14 Script de estimación por máxima verosimilitud
Al igual que con el método de los momentos, para una distribución dada, la forma de estimar por momentos es siempre la misma. Por ejemplo, para una distribución normal, los parámetros µ y σ siempre se despejan igual.
Para todas las distribuciones para las que se ha calculado el estimador de máxima verosimilitud, ya sea en clase o en los ejercicios, haz un programa en Python que tome como datos de
entrada una variable que indique el tipo de distribución (de entre aquellas para las que hemos calculado los estimadores de máxima verosimilitud, ya sea en la parte teórica o en los ejercicios),
y un array numpy que contenga los valores de la muestra y que devuelva la estimación de los parámetros correspondiente.
"""


def estimacion_maxima_verosimilitud(distribucion, valores_muestra):
    
    parametros = {}

    if distribucion == 'poisson': 
        """
            La estimación por máxima verosimilitud del parametro de Poisson es la media
        """

        return valores_muestra.mean()

    if distribucion == 'normal': 

        """
            La estimación de mu es la media muestral y el de la varianza es la varianza muestral
        """


        media_muestral  = valores_muestra.mean()
        sigma = 0
        for i in valores_muestra:
            sigma += (i - media_muestral)**2
        sigma /= len(valores_muestra)

        return {'mu': media_muestral, 'sigma': sigma}

    if distribucion == 'uniforme': #

        # En el caso de una uniforme el máximo valor de L(theta) se alcanza en el máximo valor 
        # de los xi ya que para que la derivada de la funcion de verosimilitud sea 0 
        # n debe ser 0 o theta debe ser infinito cosa que no ocurre
        
        return max(valores_muestra)

    if distribucion == 'binomial': 
        """
            La estimación de la p es la suma de todos los xi / (n * k)
        """

        return valores_muestra.sum() / (len(valores_muestra) * max(valores_muestra))
    
    if distribucion == 'geometrica':



        return len(valores_muestra) / valores_muestra.sum()

    if distribucion == 'exponencial': 

        """
            La estimación es 1/media
        """

        return 1 / valores_muestra.mean()

    if distribucion == 'rayleigh': #

        muestras_al_cuadrado = np.array([x**2 for x in valores_muestra])

        return math.sqrt(muestras_al_cuadrado.mean()/2)
           
    if distribucion == 'definida': #
        t = symbols('theta')
        n1, n2, n3, n4 = valores_muestra
        primera_derivada = n1 * (0.25/(0.5 + 0.25 * t)) + n2 * (0.25/(0.25 * t)) - n3 * (0.25/(0.25 - t * 0.25)) - n4 * (0.25/(0.25 - t * 0.25)) 
        segunda_derivada = diff(primera_derivada, t)
        soluciones_derivada = solve(primera_derivada)

        for solucion in soluciones_derivada:
            if segunda_derivada.subs(t, solucion) < 0:
                   parametros['theta'] = solucion  

        return parametros

    return None


print(estimacion_maxima_verosimilitud('definida', [1997, 32, 906, 904]))
