import sympy as sympy
from sympy import Symbol, sympify, Abs, diff
from sympy.abc import x
import numpy as np
from scipy import linalg


def biseccion(fx, Tol, Niter, a, b):

    output = {
        "type": 1,
        "method": "Bisection",
        "columns": ["iter", "a", "xm", "b", "f(xm)", "E"],
        "iterations": Niter,
        "errors": list()
    }

    datos = list()
    x = Symbol('x')
    i = 1
    cond = Tol
    error = 1.0000000

    Fun = sympify(fx)

    ea = Fun.subs(x, a)
    ea = ea.evalf()

    xm0 = 0.0
    ex_3 = 0

    xm = (a + b)/2

    ex_3 = Fun.subs(x, xm)
    ex_3 = ex_3.evalf()

    try:
        datos.append([0, '{:^15.7f}'.format(a), '{:^15.7f}'.format(
            xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(ex_3)])
        while (error > cond) and (i < Niter):
            if (ea*ex_3 < 0):
                b = xm
            else:
                a = xm

            xm0 = xm
            xm = (a+b)/2

            ex_3 = Fun.subs(x, xm)
            ex_3 = ex_3.evalf()

            error = Abs(xm-xm0)

            datos.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(
                xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(ex_3), '{:^15.7E}'.format(error)])

            i += 1
    except BaseException as e:
        if str(e) == "can't convert complex to float":
            output["errors"].append(
                "Error in data: found complex in calculations")
        else:
            output["errors"].append("Error in data: " + str(e))

        return output

    output["results"] = datos
    output["root"] = xm
    return output


def puntoFijo(X0, Tol, Niter, fx, gx):

    output = {
        "type": 1,
        "method": "Fixed point",
        "columns": ["iter", "xi", "g(xi)", "f(xi)", "E"],
        "iterations": Niter,
        "errors": list()
    }

    datos = list()
    x = sympy.Symbol('x')
    i = 1
    cond = Tol
    error = 1.000

    ex = sympify(fx)
    rx = sympify(gx)

    xP = X0
    xA = 0.0

    ea = ex.subs(x, xP)
    ea = ea.evalf()

    ra = rx.subs(x, xP)
    ra = ra.evalf()

    datos.append([0, '{:^15.7f}'.format(float(xA)), '{:^15.7f}'.format(
        float(ra)), '{:^15.7E}'.format(float(ea))])
    try:
        while((error > cond) and (i < Niter)):

            ra = rx.subs(x, xP)
            xA = ra.evalf()

            ea = ex.subs(x, xA)
            ea = ea.evalf()

            error = Abs(xA - (xP))

            xP = xA

            datos.append([i, '{:^15.7f}'.format(float(xA)), '{:^15.7f}'.format(
                float(ra)), '{:^15.7E}'.format(float(ea)), '{:^15.7E}'.format(float(error))])

            i += 1

    except BaseException as e:
        output["errors"].append("Error in data: " + str(e))
        return output

    output["results"] = datos
    output["root"] = xA
    return output


def newton(X0, Tol, Niter, fx, df):

    output = {
        "type": 1,
        "columns": ["N", "xi", "F(xi)", "E"],
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
    derf = []
    xi = X0
    f = Fun.evalf(subs={x: X0})
    derivada = DerF.evalf(subs={x: X0})
    c = 0
    Error = 100

    fn.append(f)
    xn.append(xi)
    E.append(Error)
    N.append(c)

    try:
        datos.append([c, '{:^15.7f}'.format(X0), '{:^15.7f}'.format(f)])
        while Error > Tol and f != 0 and derivada != 0 and c < Niter:
            xi = xi-f/derivada
            derivada = DerF.evalf(subs={x: xi})
            f = Fun.evalf(subs={x: xi})
            fn.append(f)
            xn.append(xi)
            c = c+1
            Error = abs(xn[c]-xn[c-1])
            derf.append(derivada)
            datos.append([c, '{:^15.7f}'.format(float(xi)), '{:^15.7E}'.format(
                float(f)), '{:^15.7E}'.format(float(Error))])
    except BaseException as e:
        if str(e) == "can't convert complex to float":
            output["errors"].append(
                "Error in data: found complex in calculations")
        else:
            output["errors"].append("Error in data: " + str(e))
        return output

    output["results"] = datos
    output["root"] = xi
    return output


def reglaFalsa(a, b, Niter, Tol, fx):

    output = {
        "type": 1,
        "method": "Regula falsi",
        "columns": ["iter", "a", "xm", "b", "f(xm)", "E"],
        "iterations": Niter,
        "errors": list()
    }

    datos = list()
    x = sympy.Symbol('x')
    i = 1
    cond = Tol
    error = 1.0000000

    ex = sympify(fx)

    xm = 0
    xm0 = 0
    ex_2 = 0
    ex_3 = 0
    ex_a = 0
    ex_b = 0

    try:
        while (error > cond) and (i < Niter):
            if i == 1:
                ex_2 = ex.subs(x, a)
                ex_2 = ex_2.evalf()
                ex_a = ex_2

                ex_2 = ex.subs(x, b)
                ex_2 = ex_2.evalf()
                ex_b = ex_2

                xm = (ex_b*a - ex_a*b)/(ex_b-ex_a)
                ex_3 = ex.subs(x, xm)
                ex_3 = ex_3.evalf()
                datos.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(
                    xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(ex_3)])
            else:

                if (ex_a*ex_3 < 0):
                    b = xm
                else:
                    a = xm

                xm0 = xm
                ex_2 = ex.subs(x, a)
                ex_2 = ex_2.evalf()
                ex_a = ex_2

                ex_2 = ex.subs(x, b)
                ex_2 = ex_2.evalf()
                ex_b = ex_2

                xm = (ex_b*a - ex_a*b)/(ex_b-ex_a)

                ex_3 = ex.subs(x, xm)
                ex_3 = ex_3.evalf()

                error = Abs(xm-xm0)
                er = sympify(error)
                error = er.evalf()
                datos.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(
                    xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(ex_3), '{:^15.7E}'.format(error)])
            i += 1
    except BaseException as e:
        if str(e) == "can't convert complex to float":
            output["errors"].append(
                "Error in data: found complex in calculations")
        else:
            output["errors"].append("Error in data: " + str(e))

        return output

    output["results"] = datos
    output["root"] = xm
    return output


def secante(fx, tol, Niter, x0, x1):
    output = {
    "type": 1,
    "method": "Secant",
    "columns": ["iter", "xi", "f(xi)", "E"],
    "errors": list()
    }

    results = list()
    x = Symbol('x')
    i = 0
    cond = tol
    error = 1.0000000

    ex = sympify(fx)

    y = x0
    ex_0 = ex
    ex_1 = ex

    try:
        while((error > cond) and (i < Niter)):
            if i == 0:
                ex_0 = ex.subs(x, x0)
                ex_0 = ex_0.evalf()
                results.append([i, '{:^15.7f}'.format(
                    float(x0)), '{:^15.7E}'.format(float(ex_0))])
            elif i == 1:
                ex_1 = ex.subs(x, x1)
                ex_1 = ex_1.evalf()
                results.append([i, '{:^15.7f}'.format(
                    float(x1)), '{:^15.7E}'.format(float(ex_1))])
            else:
                y = x1
                x1 = x1 - (ex_1*(x1 - x0)/(ex_1 - ex_0))
                x0 = y

                ex_0 = ex.subs(x, x0)
                ex_0 = ex_1.evalf()

                ex_1 = ex.subs(x, x1)
                ex_1 = ex_1.evalf()

                error = Abs(x1 - x0)

                results.append([i, '{:^15.7f}'.format(float(x1)), '{:^15.7E}'.format(
                    float(ex_1)), '{:^15.7E}'.format(float(error))])
            i += 1
    except BaseException as e:
        if str(e) == "can't convert complex to float":
            output["errors"].append(
                "Error in data: found complex in calculations")
        else:
            output["errors"].append("Error in data: " + str(e))

        return output

    output["results"] = results
    output["root"] = y
    return output


def raicesMultiples(fx, x0, tol, niter):

    output = {
        "type": 1,
        "method": "Multi Roots",
        "columns": ["iter", "xi", "f(xi)", "E"],
        "iterations": niter,
        "errors": list()
    }
    results = list()

    x = Symbol('x')

    cond = tol
    error = 1.0000

    ex = sympify(fx)

    d_ex = diff(ex, x)  # Primera derivada de Fx
    d2_ex = diff(d_ex, x)  # Segunda derivada de Fx

    xP = x0
    ex_2 = ex.subs(x, x0)
    ex_2 = ex_2.evalf()

    d_ex2 = d_ex.subs(x, x0)
    d_ex2 = d_ex2.evalf()

    d2_ex2 = d2_ex.subs(x, x0)
    d2_ex2 = d2_ex2.evalf()

    i = 0
    results.append([i, '{:^15.7E}'.format(x0), '{:^15.7E}'.format(ex_2)])
    try:
        while((error > cond) and (i < niter)):
            if(i == 0):
                ex_2 = ex.subs(x, xP)
                ex_2 = ex_2.evalf()
            else:
                d_ex2 = d_ex.subs(x, xP)
                d_ex2 = d_ex2.evalf()

                d2_ex2 = d2_ex.subs(x, xP)
                d2_ex2 = d2_ex2.evalf()

                xA = xP - (ex_2*d_ex2)/((d_ex2)**2 - ex_2*d2_ex2)

                ex_A = ex.subs(x, xA)
                ex_A = ex_A.evalf()

                error = Abs(xA - xP)
                error = error.evalf()
                er = error

                ex_2 = ex_A
                xP = xA

                results.append([i, '{:^15.7E}'.format(float(xA)), '{:^15.7E}'.format(
                    float(ex_2)), '{:^15.7E}'.format(float(er))])
            i += 1
    except BaseException as e:
        if str(e) == "can't convert complex to float":
            output["errors"].append(
                "Error in data: found complex in calculations")
        else:
            output["errors"].append("Error in data: " + str(e))

        return output

    output["results"] = results
    output["root"] = xA
    return output


def jacobi(Ma, Vb, x0, tol, niter):
    print(Ma)
    print(x0)
    print(Vb)

    output = {        
        "type": 4,
        "method": "Jacobi's Method",
        "iterations": niter,
        "errors": list(),
    }

    sX = np.size(x0)
    xA = np.zeros((sX,1))

    A = np.matrix(Ma)

    b = np.array(Vb)
    s = b.size
    b = np.reshape(b,(s,1))

    D = np.diag(np.diag(A))
    L = -1*np.tril(A)+D
    U = -1*np.triu(A)+D
    LU = L+U

    T = np.linalg.inv(D) @ LU
    C = np.linalg.inv(D) @ b


    output["t"] = T 
    output["c"]= C

    print(T)
    print(C)
    return output

def gaussSeidel(Ma, Vb, x0, tol, niter):
    output = {
        "type": 4,
        "method": "Gauss-Seidel's Method",
        "iterations": niter
    }
    
    sX = np.size(x0)
    xA = np.zeros((sX,1))
    
    A = np.matrix(Ma)
    
    b = np.array(Vb)
    s = b.size
    b = np.reshape(b,(s,1))

    D = np.diag(np.diag(A))
    L = -1*np.tril(A)+D
    U = -1*np.triu(A)+D
    
    T = np.linalg.inv(D-L) @ U
    C = np.linalg.inv(D-L) @ b
    
    xP = x0
    E = 1000
    cont = 0

    steps = {'Step 0': np.copy(xA)}
    while(E > tol and cont < niter  ):
        xA = T@xP + C
        E = np.linalg.norm(xP - xA)
        xP = xA
        cont = cont + 1
        steps[f'Step {cont+1}'] = np.copy(xA)

    niter = cont
    error = E

    print("T", T)
    print("C", C)
    print("steps", steps)

def splineLineal(X,Y):
    output = {
        "type": 8,
        "method": "Linear Tracers",
        "errors" : list()
    }
    X = np.array(X)
    Y = np.array(Y)
    n = X.size
    m = 2*(n-1)
    A = np.zeros((m,m))
    b = np.zeros((m,1))
    Coef = np.zeros((n-1,2))
    i = 0
    #Interpolating condition
    try:
        while i < X.size-1:
            A[i+1,[2*i+1-1,2*i+1]]= [X[i+1],1] 
            b[i+1]=Y[i+1]
            i = i+1

        A[0,[0,1]] = [X[0],1] 
        b[0] = Y[0]
        i = 1
        #Condition of continuity
        while i < X.size-1:
            A[X.size-1+i,2*i-2:2*i+2] = np.hstack((X[i],1,-X[i],-1))
            b[X.size-1+i] = 0
            i = i+1

        Saux = linalg.solve(A,b)
        #Order Coefficients
        i = 0
        while i < X.size-1:
            Coef[i,:] = [Saux[2*i],Saux[2*i+1]]
            i = i+1
        
    except BaseException as e:  
        output["errors"].append("Error in data: " + str(e))
        return output
    
    output["results"] = Coef
    return output

# funcional
def sor(Ma, Vb, x0, w, tol, niter):
    output = {
        "type": 4,
        "method": "SOR(Relaxation) Method",
        "iterations": niter
    }
    
    sX = np.size(x0)
    xA = np.zeros((sX,1))
    
    A = np.matrix(Ma)
    
    b = np.array(Vb)
    s = b.size
    b = np.reshape(b,(s,1))

    D = np.diag(np.diag(A))
    L = -1*np.tril(A)+D
    U = -1*np.triu(A)+D
    
    T = np.linalg.inv(D-(w*L)) @ (((1-w)*D)+(w*U))
    C = (w*np.linalg.inv(D-(w*L))) @ b
    
    xP = x0
    E = 1000
    cont = 0
    steps = {'Step 0': np.copy(xA)}
    while(E > tol and cont < niter):
        xA = T@xP + C
        E = np.linalg.norm(xP - xA)
        xP = xA
        cont = cont + 1
        steps[f'Step {cont+1}'] = np.copy(xA)

    nIter = cont
    error = E
    
    output["results"] = steps
    output["E"] = error
    output["Iteration"] = nIter
    
    print(T)
    print(C)
    print("steps", steps)
    
    return output

