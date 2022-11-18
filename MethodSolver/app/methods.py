import sympy as sympy


from sympy.abc import x
import pandas as pd
import math
import numpy as np
import numpy.matlib
from numpy.lib import scimath
from sympy import sympify, diff, evalf, Abs, Symbol, Subs
from scipy import linalg

def biseccion(fx, Tol, Niter, Xs, Xi):

    output = {
        "type": 1,
        "columns": ["N","xm","F(xm)","E" ],
        "errors": list()
    }

    datos = list()
    x = sympy.Symbol('x')
    Fun = sympy.sympify(fx, convert_xor=True)


    fm = []
    E = []
    xAux = Xi
    fi = Fun.evalf(subs={x: Xi})
    xAux = Xs
    fs = Fun.evalf(subs={x: Xs})
    iter = []
    Xms = []

    if fi == 0:
        s = Xi
        E = 0
        print(Xi, "es raiz de f(x)")
    elif fs == 0:
        s = Xs
        E = 0
        print(Xs, "es raiz de f(x)")
    elif fs*fi < 0:
        c = 0
        Xm = (Xi+Xs)/2
        xAux = Xm
        fe = Fun.evalf(subs={x: xAux})
        fm.append(fe)
        E.append(100)
        iter.append(c)
        Xms.append(xAux)

        while E[c] > Tol and fe != 0 and c < Niter:
            if fi*fe < 0:
                Xs = Xm
                xAux = Xs
                fs = Fun.evalf(subs={x: xAux})
            else:
                Xi = Xm
                xAux = Xi
                fs = Fun.evalf(subs={x: xAux})
            Xa = Xm
            Xm = (Xi+Xs)/2
            xAux = Xm
            fe = Fun.evalf(subs={x: xAux})
            fm.append(fe)
            Error = abs(Xm-Xa)
            E.append(Error)
            c = c+1
            iter.append(c)
            Xms.append(xAux)

            # Se guardan los datos            
        
        datos.append([iter, Xms, fm, E])

        if fe == 0:
            s = xAux
            print(s, "es raiz de f(x)")
        elif Error < Tol:
            s = xAux
            #print("//", "Fm:", fm, "//")
            #print("//", "Error:", fm, "//")
            #print(s, "es una aproximacion de un raiz de f(x) con una tolerancia", Tol)
            output["results"] = datos
            output["root"] = s
            return output

        else:
            s = xAux
            print("Fracaso en ", Niter, " iteraciones ")
            return 

    else:
        print("El intervalo es inadecuado")


def puntoFijo(X0, Tol, Niter, fx, gx):

    output = {
        "type": 1,
        "columns": ["N","xi","F(xi)", "G(xi)", "E" ],
        "errors": list()
    }

    datos = list()
    x = sympy.Symbol('x')
    Fun = sympy.sympify(fx, convert_xor=True)
    GFun = sympy.sympify(gx, convert_xor=True) 

    gn = []
    fn = []
    xn = []
    E = []
    N = []
    xi = X0
    f = Fun.evalf(subs={x: X0})
    c = 0
    Error = 100

    fn.append(f)
    xn.append(xi)
    E.append(Error)
    N.append(c)
    while Error > Tol and f != 0 and c < Niter:
        xi = GFun.evalf(subs={x: xi})
        fe = Fun.evalf(subs={x: xi})
        gxFun = GFun.evalf(subs={x: xi})
        gn.append(gxFun)
        fn.append(fe)
        xn.append(xi)
        c = c+1
        Error = abs(xn[c]-xn[c-1])
        N.append(c)
        E.append(Error)
    datos.append([N, xn, fn, gn, E])
    if fe == 0:
        s = xi
        print(s, "es raiz de f(x)")
    elif Error < Tol:
        s = xi
        print(s, "es una aproximacion de un raiz de f(x) con una tolerancia", Tol)
        print("N", N)
        print("xn", xn)
        print("fn", fn)
        print("Error", E)
        output["results"] = datos
        output["root"] = s
        return output
    else:
        s = xi
        print("Fracaso en ", Niter, " iteraciones ")


def newton(X0, Tol, Niter, Fun, df):
    fn = []
    xn = []
    E = []
    N = []
    x = X0
    f = eval(Fun)
    derivada = eval(df)
    c = 0
    Error = 100
    fn.append(f)
    xn.append(x)
    E.append(Error)
    N.append(c)
    while Error > Tol and f != 0 and derivada != 0 and c < Niter:
        x = x-f/derivada
        derivada = eval(df)
        f = eval(Fun)
        fn.append(f)
        xn.append(x)
        c = c+1
        Error = abs(xn[c]-xn[c-1])
        N.append(c)
        E.append(Error)
        if f == 0:
            s = x
            print(s, "es raiz de f(x)")
        elif Error < Tol:
            s = x
            print(s, "es una aproximacion de un raiz de f(x) con una tolerancia", Tol)
            print("N", N)
            print("xn", xn)
            print("fn", fn)
            print("Error", E)
        else:
            s = x
            print("Fracaso en ", Niter, " iteraciones ")


def reglaFalsa(Xi, Xf, Niter, Tol, Fun):
    # iniciar variables
    solucion = None
    contador = 0
    error = 101
    # Evaluar si la raiz esta dentro del intervalo
    if Fun(Xi)*Fun(Xf) <= 0:
        while contador <= Niter and error >= Tol:
            contador += 1
            solucion = Xf-((Fun(Xf)*(Xf-Xi))/(Fun(Xf)-Fun(Xi)))
            error = abs((solucion-Xi)/solucion)*100
            # redeefinir el nuevo intervalo
            if Fun(Xi)*Fun(solucion) >= 0:
                Xi = solucion
            else:
                Xf = solucion
        print(solucion)
        print(contador)
        print(error)
    else:
        print("no hay solucion")
