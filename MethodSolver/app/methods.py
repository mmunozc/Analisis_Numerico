import sympy as sympy


from sympy.abc import x
import pandas as pd
import numpy as np
import math


def biseccion(fx, tol, niter, xs, xi):
    x = sympy.Symbol('x')

    Fun = sympy.sympify(fx, convert_xor=True)
    print
    Tol = float(tol)
    Niter = int(niter)
    Xs = float(xs)
    Xi = float(xi)
    
    fm = []
    E = []
    xAux = Xi
    fi = Fun.evalf(subs={x: Xi})
    xAux = Xs
    fs = Fun.evalf(subs={x: Xs})
    
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
        if fe == 0:
            s = xAux
            print(s, "es raiz de f(x)")
        elif Error < Tol:
            s = xAux
            #print("//", "Fm:", fm, "//")
            #print("//", "Error:", fm, "//")
            #print(s, "es una aproximacion de un raiz de f(x) con una tolerancia", Tol)

            return fm, Error
        else:
            s = xAux
            print("Fracaso en ", Niter, " iteraciones ")
    else:
        print("El intervalo es inadecuado")




