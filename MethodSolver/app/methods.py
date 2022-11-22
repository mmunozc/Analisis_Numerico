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



#funcional
def pivParcial(a,b):
    datos=list()
    ab=np.concatenate((a,b), axis=1)
    #pivoteo por fila parcial
    tamano=np.shape(ab)
    n=tamano[0]
    m=tamano[1]

    for i in range(0, n-1, 1):
        columna=abs(ab[i:,i])
        dondemax=np.argmax(columna)
        #intercambio de fila
        if dondemax!=0:
            temporal=np.copy(ab[i,:])
            ab[i,:]=ab[dondemax+i, :]
            ab[dondemax,:]=temporal
    datos.append(ab)
    return datos



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
        output["root"] = s
        return output
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


def newton(X0, Tol, Niter, fx, df):

    output = {
        "type": 1,
        "columns": ["N","xi","F(xi)", "E" ],
        "errors": list()
    }

    datos = list()
    x = sympy.Symbol('x')
    Fun = sympy.sympify(fx, convert_xor=True)
    DerF = sympy.sympify(df, convert_xor=True) 


    fn = []
    xn = []
    E = []
    N = []
    derf=[]
    xi = X0
    f = Fun.evalf(subs={x: X0})
    derivada = DerF.evalf(subs={x: X0})
    c = 0
    Error = 100

    fn.append(f)
    xn.append(xi)
    E.append(Error)
    N.append(c)

    while Error > Tol and f != 0 and derivada != 0 and c < Niter:
        xi = xi-f/derivada
        derivada = DerF.evalf(subs={x: xi})
        f = Fun.evalf(subs={x: xi})
        fn.append(f)
        xn.append(xi)
        c = c+1
        Error = abs(xn[c]-xn[c-1])
        derf.append(derivada)
        N.append(c)
        E.append(Error)
        
    datos.append([N, xn, fn, E])
    
    if f == 0:
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


#Arreglar da 101 iteraciones
def reglaFalsa(Xi, Xf, Niter, Tol, fx):

    output = {
        "type": 1,
        "columns": ["N","xm","F(xm)", "E" ],
        "errors": list()
    }

    datos = list()
    x = sympy.Symbol('x')
    Fun = sympy.sympify(fx, convert_xor=True)

    fn = []
    xm = []
    E = []
    N = []

    fxi = Fun.evalf(subs={x: Xi})
    fxf = Fun.evalf(subs={x: Xf})

    # iniciar variables
    solucion = None
    c = 0
    error = 1.000000

    fn.append(fxi)
    xm.append(Xi)
    E.append(error)
    N.append(c)

    # Evaluar si la raiz esta dentro del intervalo
    if fxi*fxf <= 0:
        while c <= Niter and error >= Tol:
            c += 1
            solucion = (Xf-(fxf*(Xf-Xi))/(fxf-fxi))
            error = abs((solucion-Xi)/solucion)*100
            # redeefinir el nuevo intervalo
            fxs = Fun.evalf(subs={x: solucion})
            if fxi*fxs >= 0:
                Xi = solucion
            else:
                Xf = solucion
        
            fn.append(solucion)
            xm.append(Xi)
            E.append(error)
            N.append(c)
        print(solucion)
        print(N)
        print(E)
    else:
        print("no hay solucion")


def secante(f, x1, x2, tol):
    error = 1e3
    n = 0
    x3 = 0
    while error > tol:
        x3 = x-((x2-x1)/(f(x2)-f(x1)))*f(x1)
        x1=x2
        x2=x3
        error = abs(f(x3))
        n+=1
    pass


def raphson(f, fp ,xi ,tol ,Niter):

    datos = list()
    x = sympy.Symbol('x')
    Fun = sympy.sympify(f, convert_xor=True)
    FunP = sympy.sympify(fp, convert_xor=True)
    
    fx = Fun.evalf(subs={x: xi})
    fxp = FunP.evalf(subs={x: xi})

    fxi = []
    Xi = []
    E = []
    N = []
    
    for k in range(Niter):
        xold=xi
        xi=xi-fx/fxp
        e=abs((xi-xold)/xi)
        if e < tol:
            break
        #datos.append(k, xi, fx, e)
        fxi.append(xold)
        Xi.append(xi)
        E.append(e)
        N.append(k)
    
    print(E)

#funcional
def luDirecta(n, a):
    matriz=a
    datos=list()
    u=np.zeros([n,n])
    l=np.zeros([n,n])

    for r in range(0, n):
        for c in range(0, n):
            u[r][c]=matriz[r][c]
    
#operacion para hacer 0 debajo de la diagonal
    for k in range(0,n):
        for r in range(0,n):
            if k==r:
                l[k,r]=1
            if k<r:
                factor=(matriz[r][k]/matriz[k][k])
                l[r][k]=factor
                for c in range(0,n):
                    matriz[r][c]=matriz[r][c]-(factor*matriz[k][c])
                    u[r][c]=matriz[r][c]
    datos.append([l, u])
    #l=np.transpose(l)
    #u=np.transpose(u)
    return datos
            
