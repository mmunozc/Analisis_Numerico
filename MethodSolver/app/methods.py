import sympy as sympy
from sympy import Symbol, sympify, Abs, diff
from sympy.abc import x
import numpy as np
from scipy import linalg

#fx/Fun = función
#x0= valor incial
#a,b = puntos intervalo
#Tol = toleracia
#f(letra/s)=función evaluada en esa variable


def biseccion(fx, Tol, Niter, a, b):
    output = {
        "columns": ["iter", "a", "xm", "b", "f(xm)", "E"],
        "iterations": Niter,
        "errors": list()
    }

    # Configuraciones iniciales
    datos = list()
    x = Symbol('x')
    i = 1
    error = 1.0000000
    Fun = sympify(fx)

    Fa = Fun.subs(x, a) # Funcion evaluada en a
    Fa = Fa.evalf()

    xm0 = 0.0
    Fxm = 0

    xm = (a + b)/2 # Punto intermedio

    Fxm = Fun.subs(x, xm) # Funcion evaluada en Xm
    Fxm = Fxm.evalf()

    try:
        datos.append([0, '{:^15.7f}'.format(a), '{:^15.7f}'.format(xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(Fxm)]) # Datos con formato dado
        while (error > Tol) and (i < Niter): # Se repite hasta que el intervalo sea lo pequeño que se desee
            if (Fa*Fxm < 0): # Se elecciona un intervalo inicial, donde el valor de la funcion cambie de signo en [a,b]
                b = xm
            else:
                a = xm # Cambia de signo en [m,b]

            xm0 = xm
            xm = (a+b)/2 # Se calcula el punto intermedio del intervalo - Divide el intervalo a la mitadd

            Fxm = Fun.subs(x, xm)
            Fxm = Fxm.evalf() # Se evalua el punto intermedio en la funcion

            error = Abs(xm-xm0) # Se calcula el error

            datos.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(xm), 
                            '{:^15.7f}'.format(b), '{:^15.7E}'.format(Fxm), '{:^15.7E}'.format(error)]) # Se van agregando las soluciones con el formato deseado

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

#Completar los comentarios y explicar mas a fondo
def puntoFijo(X0, Tol, Niter, fx, gx):

    output = {
        "columns": ["iter", "xi", "g(xi)", "f(xi)", "E"],
        "iterations": Niter,
        "errors": list()
    }

    #configuración inicial
    datos = list()
    x = sympy.Symbol('x')
    i = 1
    Tol
    error = 1.000

    Fx = sympify(fx)
    Gx = sympify(gx)

    # Iteracion 0
    xP = X0 # Valor inicial (Punto evaluacion)
    xA = 0.0

    Fa = Fx.subs(x, xP) # Funcion evaluada en el valor inicial
    Fa = Fa.evalf()

    Ga = Gx.subs(x, xP) # Funcion G evaluada en el valor inicial
    Ga = Ga.evalf()

    datos.append([0, '{:^15.7f}'.format(float(xA)), '{:^15.7f}'.format(
        float(Ga)), '{:^15.7E}'.format(float(Fa))])
    try:
        while((error > Tol) and (i < Niter)): # Se repite hasta que el error sea menor a la tolerancia

            # Se evalua el valor inicial en G, para posteriormente evaluar este valor en la funcion F siendo-> Xn=G(x) y F(xn) = F(G(x))
            Ga = Gx.subs(x, xP) # Funcion G evaluada en el punto de inicial
            xA = Ga.evalf()

            Fa = Fx.subs(x, xA)# Funcion evaluada en el valor de la evalucacion de G
            Fa = Fa.evalf()

            error = Abs(xA - (xP)) # Se calcula el error 

            xP = xA # Nuevo punto de evaluacion (Punto inicial)

            datos.append([i, '{:^15.7f}'.format(float(xA)), '{:^15.7f}'.format(
                float(Ga)), '{:^15.7E}'.format(float(Fa)), '{:^15.7E}'.format(float(error))])

            i += 1

    except BaseException as e:
        output["errors"].append("Error in data: " + str(e))
        return output

    output["results"] = datos
    output["root"] = xA
    return output

#Comentado
def newton(x0, Tol, Niter, fx, df):

    output = {
        "columns": ["N", "xi", "F(xi)", "E"],
        "errors": list()
    }

    #configuración inicial
    datos = list()
    x = sympy.Symbol('x')
    Fun = sympify(fx)
    DerF = sympify(df)


    xn = []
    derf = []
    xi = x0 # Punto de inicio
    f = Fun.evalf(subs={x: x0}) #función evaluada en x0
    derivada = DerF.evalf(subs={x: x0}) #función derivada evaluada en x0
    c = 0
    Error = 100
    xn.append(xi)

    try:
        datos.append([c, '{:^15.7f}'.format(x0), '{:^15.7f}'.format(f)])

        # Al evaluar la derivada en el punto inicial,se busca que sea diferente de 9, ya que al serlo nos encontramos en un punto de inflexion
        #(No se puede continuar ya que la tangente es horinzontal)
        while Error > Tol and f != 0 and derivada != 0 and c < Niter: # El algoritmo converge o se alcanzo limite de iteraciones fijado

            xi = xi-f/derivada # Estimacion del siguiente punto aproximado a la raiz (nuevo valor inicial)
            derivada = DerF.evalf(subs={x: xi}) # Evaluacion de la derivada con el nuevo valor inicial (xi)
            f = Fun.evalf(subs={x: xi}) # Evaluacion de la derivada con el nuevo valor inicial (xi)
            xn.append(xi)
            c = c+1
            Error = abs(xn[c]-xn[c-1]) # Se reduce entre cada iteracion (Representado por el tramo)
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
        "columns": ["iter", "a", "xm", "b", "f(xm)", "E"],
        "iterations": Niter,
        "errors": list()
    }

    #configuración inicial
    datos = list()
    x = sympy.Symbol('x')
    i = 1
    cond = Tol
    error = 1.0000000

    Fun = sympify(fx)

    xm = 0
    xm0 = 0
    Fx_2 = 0
    Fx_3 = 0
    Fa = 0
    Fb = 0

    try:
        while (error > cond) and (i < Niter):
            if i == 1:
                Fx_2 = Fun.subs(x, a)
                Fx_2 = Fx_2.evalf()
                Fa = Fx_2

                Fx_2 = Fun.subs(x, b)
                Fx_2 = Fx_2.evalf()
                Fb = Fx_2

                xm = (Fb*a - Fa*b)/(Fb-Fa)
                Fx_3 = Fun.subs(x, xm)
                Fx_3 = Fx_3.evalf()
                datos.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(
                    xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(Fx_3)])
            else:

                if (Fa*Fx_3 < 0):
                    b = xm
                else:
                    a = xm

                xm0 = xm
                Fx_2 = Fun.subs(x, a)
                Fx_2 = Fx_2.evalf()
                Fa = Fx_2

                Fx_2 = Fun.subs(x, b)
                Fx_2 = Fx_2.evalf()
                Fb = Fx_2

                xm = (Fb*a - Fa*b)/(Fb-Fa)

                Fx_3 = Fun.subs(x, xm)
                Fx_3 = Fx_3.evalf()

                error = Abs(xm-xm0)
                er = sympify(error)
                error = er.evalf()
                datos.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(
                    xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(Fx_3), '{:^15.7E}'.format(error)])
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

#Comentado
def secante(fx, tol, Niter, x0, x1):
    output = {
        "columns": ["iter", "xi", "f(xi)", "E"],
        "errors": list()
    }

    results = list()
    x = Symbol('x')
    i = 0
    cond = tol
    error = 1.0000000

    Fun = sympify(fx)

    y = x0
    Fx0 = Fun
    Fx1 = Fun

    try:
        while((error > cond) and (i < Niter)):
            if i == 0:
                Fx0 = Fun.subs(x, x0) #Evaluacion en el valor inicial X0
                Fx0 = Fx0.evalf()
                results.append([i, '{:^15.7f}'.format(float(x0)), '{:^15.7E}'.format(float(Fx0))])
            elif i == 1:
                Fx1 = Fun.subs(x, x1)#Evaluacion en el valor inicial X1
                Fx1 = Fx1.evalf()
                results.append([i, '{:^15.7f}'.format(float(x1)), '{:^15.7E}'.format(float(Fx1))])
            else:
                y = x1 
                # Se calcula la secante
                x1 = x1 - (Fx1*(x1 - x0)/(Fx1 - Fx0)) # Punto de corte del intervalo usando la raiz de la secante, (xi+1)
                x0 = y

                Fx0 = Fun.subs(x, x0) #Evaluacion en el valor inicial X0
                Fx0 = Fx1.evalf() 

                Fx1 = Fun.subs(x, x1)#Evaluacion en el valor inicial X1
                Fx1 = Fx1.evalf()

                error = Abs(x1 - x0) # Tramo

                results.append([i, '{:^15.7f}'.format(float(x1)), '{:^15.7E}'.format(float(Fx1)), '{:^15.7E}'.format(float(error))])
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

#Comentado
def raicesMultiples(fx, x0, tol, niter):

    output = {
        "columns": ["iter", "xi", "f(xi)", "E"],
        "iterations": niter,
        "errors": list()
    }
    
    # Configuraciones inciales
    results = list()
    x = Symbol('x')
    cond = tol
    error = 1.0000
    ex = sympify(fx)

    d_ex = diff(ex, x)  # Primera derivada de Fx
    d2_ex = diff(d_ex, x)  # Segunda derivada de Fx


    xP = x0
    ex_2 = ex.subs(x, x0)  # Funcion evaluada en x0
    ex_2 = ex_2.evalf() 

    d_ex2 = d_ex.subs(x, x0) # primera derivada evaluada en x0
    d_ex2 = d_ex2.evalf()

    d2_ex2 = d2_ex.subs(x, x0) # segunda derivada evaluada en x0
    d2_ex2 = d2_ex2.evalf()

    i = 0
    results.append([i, '{:^15.7E}'.format(x0), '{:^15.7E}'.format(ex_2)]) # Datos con formato dado
    try:
        while((error > cond) and (i < niter)): # Se repite hasta que el intervalo sea lo pequeño que se desee
            if(i == 0):
                ex_2 = ex.subs(x, xP) # Funcion evaluada en xp
                ex_2 = ex_2.evalf()
            else:
                d_ex2 = d_ex.subs(x, xP) # Funcion evaluada en xp
                d_ex2 = d_ex2.evalf()

                d2_ex2 = d2_ex.subs(x, xP) # Funcion evaluada en xp
                d2_ex2 = d2_ex2.evalf()

                xA = xP - (ex_2*d_ex2)/((d_ex2)**2 - ex_2*d2_ex2) # Método de Newton-Raphson modificado

                ex_A = ex.subs(x, xA) # Funcion evaluada en xA
                ex_A = ex_A.evalf()

                error = Abs(xA - xP)
                error = error.evalf() # Se calcula el error
                er = error

                ex_2 = ex_A #se establece la nueva aproximación
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

#Metodos iterativos

#Comentado
def jacobi(Ma, Vb, x0, tol, niter):
    output = {
        "iterations": niter,
        "errors": list(),
    }

    sX = np.size(x0)

    A = np.matrix(Ma)

    b = np.array(Vb)
    s = b.size
    b = np.reshape(b, (s, 1)) #Rehace el tamaño del vector b

    D = np.diag(np.diag(A)) #saca la diagonal de la matriz A
    L = -1*np.tril(A)+D #saca la matriz Lower de la matriz A
    U = -1*np.triu(A)+D #Saca la matriz Upper de la matriz A
    LU = L+U

    T = np.linalg.inv(D) @ LU #Obtiene la matriz de Transicion multiplicando el inverso de D por la matriz LU
    C = np.linalg.inv(D) @ b #Obtiene la matriz de coeficientes multiplicando el inverso de la matriz de D por la matriz b

    output["t"] = T
    output["c"] = C

    resultado={"t":T,
                "c":C}
    return resultado

#Comentado
def gaussSeidel(Ma, Vb, x0, tol, niter):
    iteraciones=[]
    informacion=[]
    error=[]

    sX = np.size(x0)
    xA = np.zeros((sX, 1))

    A = np.matrix(Ma)

    b = np.array(Vb)
    s = b.size
    b = np.reshape(b, (s, 1)) #Rehace el tamaño del vector b

    D = np.diag(np.diag(A)) #saca la diagonal de la matriz A
    L = -1*np.tril(A)+D #saca la matriz Lower de la matriz A
    U = -1*np.triu(A)+D #Saca la matriz Upper de la matriz A

    T = np.linalg.inv(D-L) @ U #Obtiene la matriz de Transicion multiplicando el inverso de D-L por la matriz U
    tFinal=max(abs(np.linalg.eigvals(T)))
    C = np.linalg.inv(D-L) @ b #Obtiene la matriz Coeficientes multiplicando el inverso de D-L por la matriz b

    xP = x0
    E = 1000
    cont = 0

    steps = {'Step 0': np.copy(xA)}
    while(E > tol and cont < niter):
        xA = T@xP + C #Multiplica la matriz de transicion con la matriz inicial mas la de coeficientes
        E = np.linalg.norm(xP - xA) #Obtiene el error sacandole la norma a la resta de la matriz inicial con la matriz A    
        xP = xA
        cont = cont + 1
        iteraciones.append(cont)
        informacion.append(xA)
        error.append(E)
        steps[f'Step {cont+1}'] = np.copy(xA)

    datos=zip(iteraciones, error, informacion)
    resultado={"t":T,
                "c":C,
                "esp":tFinal,
                "informacion":datos}
    return resultado

#Comentado
def sor(Ma, Vb, x0, w, tol, niter):

    A = np.matrix(Ma)

    b = np.array(Vb)
    s = b.size
    b = np.reshape(b, (s, 1)) #Rehace el tamaño del vector b

    D = np.diag(np.diag(A)) #saca la diagonal de la matriz A
    L = -1*np.tril(A)+D #saca la matriz Lower de la matriz A
    U = -1*np.triu(A)+D #Saca la matriz Upper de la matriz A

    T = np.linalg.inv(D-L) @ U #Obtiene la matriz de Transicion multiplicando el inverso de D-L por la matriz U
    tFinal=max(abs(np.linalg.eigvals(T)))
    C = np.linalg.inv(D-L) @ b #Obtiene la matriz Coeficientes multiplicando el inverso de D-L por la matriz b

    iteraciones=[]
    informacion=[]
    cumple=False
    n=len(Ma)
    k=0

    while(not cumple and k<niter):
        xk1=np.zeros(n)
        for i in range(n):
            s1=np.dot(Ma[i][:i],xk1[:i]) #Multiplica los valores de la Matriz A hasta el final de la matriz xk1
            s2=np.dot(Ma[i][i+1:], x0[i+1:])# Multiplica la matrizA con el vector de inicio
            xk1[i]=(Vb[i]-s1-s2)/Ma[i][i]*w+(1-w)*x0[i] #Hace las operaciones para obtener el resultado del metodo
        norma=np.linalg.norm(x0-xk1)
        x0=xk1 #actualiza los valores para el proximo ciclo
        print('Iteracion:{}->{} norma {}'.format(k, xk1, norma))
        iteraciones.append(k)
        informacion.append(xk1)
        cumple=norma<tol
        k+=1
    
    if k<niter:
        datos=zip(iteraciones, informacion) #guarda el contador, informacion
        resultado={"solucion":x0,
                    "t":T,
                    "c":C,
                    "esp":tFinal,
                    "informacion":datos}
        return resultado
    else:
        return "el sistem no converge"

#Metodos interpolacion

def splineLineal(X, Y):
    output = {
        "errors": list()
    }
    X = np.array(X)
    Y = np.array(Y)
    n = X.size
    m = 2*(n-1) #factor m para la formula de los trazadores
    A = np.zeros((m, m))
    b = np.zeros((m, 1))
    Coef = np.zeros((n-1, 2))
    i = 0
    # Interpolating condition
    try:
        while i < X.size-1:
            A[i+1, [2*i+1-1, 2*i+1]] = [X[i+1], 1]  #Realiza la formula de interpolacion de superficie 
            b[i+1] = Y[i+1]
            i = i+1

        A[0, [0, 1]] = [X[0], 1]
        b[0] = Y[0]
        i = 1
        # Condition of continuity
        while i < X.size-1:
            A[X.size-1+i, 2*i-2:2*i+2] = np.hstack((X[i], 1, -X[i], -1))
            b[X.size-1+i] = 0
            i = i+1

        Saux = linalg.solve(A, b)
        # Order Coefficients
        i = 0
        while i < X.size-1:
            Coef[i, :] = [Saux[2*i], Saux[2*i+1]]
            i = i+1

    except BaseException as e:
        output["errors"].append("Error in data: " + str(e))
        return output

    output["results"] = Coef
    return output

def splineCuadratica(X, Y):
    output = {
        "errors": list()
    }
    X = np.array(X)
    Y = np.array(Y)
    n = X.size
    m = 3*(n-1) # Factor m para la formula de los trazadores, como es cuadratica se multiplica por 3
    A = np.zeros((m, m))
    b = np.zeros((m, 1))
    Coef = np.zeros((n-1, 3))
    i = 0
    try:
        # Interpolating condition
        while i < X.size-1:

            A[i+1, 3*i:3*i+3] = np.hstack((X[i+1]**2, X[i+1], 1)) #Realiza la formula de interpolacion de superficie 
            b[i+1] = Y[i+1]
            i = i+1

        A[0, 0:3] = np.hstack((X[0]**2, X[0]**1, 1))
        b[0] = Y[0]
        # Condition of continuity
        i = 1
        while i < X.size-1:
            A[X.size-1+i, 3*i-3:3*i +
                3] = np.hstack((X[i]**2, X[i], 1, -X[i]**2, -X[i], -1))
            b[X.size-1+i] = 0
            i = i+1
        # Condition of smoothness
        i = 1
        while i < X.size-1:
            A[2*n-3+i, 3*i-3:3*i+3] = np.hstack((2*X[i], 1, 0, -2*X[i], -1, 0))
            b[2*n-3+i] = 0
            i = i+1
        A[m-1, 0] = 2
        b[m-1] = 0

        Saux = linalg.solve(A, b)
        # Order Coefficients
        i = 0
        j = 0
        while i < n-1:
            Coef[i, :] = np.hstack((Saux[j], Saux[j+1], Saux[j+2]))
            i = i+1
            j = j + 3
    except BaseException as e:
        output["errors"].append("Error in data: " + str(e))
        return output

    output["results"] = Coef
    return output

def splineCubica(X, Y):
    output = {
        "errors": list()
    }
    X = np.array(X)
    Y = np.array(Y)
    n = X.size
    m = 4*(n-1)
    A = np.zeros((m,m))
    b = np.zeros((m,1))
    Coef = np.zeros((n-1,4))
    i = 0
    try:
        #Interpolating condition
        while i < X.size-1:
            
            A[i+1,4*i:4*i+4]= np.hstack((X[i+1]**3,X[i+1]**2,X[i+1],1)) 
            b[i+1]=Y[i+1]
            i = i+1

        A[0,0:4] = np.hstack((X[0]**3,X[0]**2,X[0]**1,1))
        b[0] = Y[0]
        #Condition of continuity
        i = 1
        while i < X.size-1:
            A[X.size-1+i,4*i-4:4*i+4] = np.hstack((X[i]**3,X[i]**2,X[i],1,-X[i]**3,-X[i]**2,-X[i],-1))
            b[X.size-1+i] = 0
            i = i+1
        #Condition of smoothness
        i = 1
        while i < X.size-1:
            A[2*n-3+i,4*i-4:4*i+4] = np.hstack((3*X[i]**2,2*X[i],1,0,-3*X[i]**2,-2*X[i],-1,0))
            b[2*n-3+i] = 0
            i = i+1
        
        #Concavity condition
        i = 1
        while i < X.size-1:
            A[3*n-5+i,4*i-4:4*i+4] = np.hstack((6*X[i],2,0,0,-6*X[i],-2,0,0))
            b[n+5+i] = 0
            i = i+1
        
        #Boundary conditions  
        A[m-2,0:2]=[6*X[0],2];
        b[m-2]=0;
        A[m-1,m-4:m-2]=[6*X[X.size-1],2];
        b[m-1]=0;
        
        Saux = linalg.solve(A,b)
        #Order Coefficients
        i = 0
        j = 0
        while i < n-1:
            Coef[i,:] = np.hstack((Saux[j],Saux[j+1],Saux[j+2],Saux[j+3]))
            i = i+1
            j = j + 4
    except BaseException as e:  
        output["errors"].append("Error in data: " + str(e))
        return output

    output["results"] = Coef
    return output

def vandermonde(a,b):
    copiaB=np.copy(b)
    longitudMatriz=len(a)
    matrizVandermonde=np.vander(a) #Obtiene la matriz vandermonde con la matrizA
    coeficientes=np.linalg.solve(matrizVandermonde, copiaB) #Encuentra la Matriz A con vector B
            
    print(coeficientes)
    x=sympy.Symbol('x')
    polinomio=0
    for i in range(0, longitudMatriz, 1): #ciclo para asignarle las x y la potencias al polinomio
        potencia=(longitudMatriz-1)-i
        termino=coeficientes[i]*(x**potencia)
        polinomio=polinomio+termino

    print(polinomio)
    datos={
        "matriz":matrizVandermonde,
        "coeficientes":coeficientes,
        "polinomio":polinomio,
    }

    return datos

#Corregir
def newtonDivDif(X, Y):
    output = {}

    X = np.array(X)
    n = X.size

    Y = np.array(Y)

    D = np.zeros((n,n))

    D[:,0]=Y.T
    for i in range(1,n):
        aux0 = D[i-1:n,i-1]
        aux = np.diff(aux0)
        aux2 = X[i:n] - X[0:n-i]
        D[i:n,i] = aux/aux2.T  
        Coef = np.diag(D)
    output["D"] = D
    output["Coef"] = Coef
    print(D)
    print(Coef)

    return output