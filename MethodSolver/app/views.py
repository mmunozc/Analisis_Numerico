from django.http import HttpResponseRedirect
from django.shortcuts import render
import numpy as np
from app.methods import *
from app.methods import *
from app.methods import *
import numpy.matlib
# Create your views here.


def infoView(request):
    return render(request, 'home.html')

def menuView(request):
    return render(request, 'metodosMenu.html')

def iterativeMethodsView(request):
    return render(request, './metodosPage/metodos-iterativos.html')
def splineMethodsView(request):
    return render(request, './metodosPage/metodos-spline.html')

def graficaView(request):
    return render(request, 'grafica.html')

def luDirectaView(request):
    datos=()
    matriz=[]
    if request.method=='POST':
        tamañoMatriz=request.POST["tamañoArreglo"]
        matriz=np.zeros((tamañoMatriz, tamañoMatriz))
        for i in range(tamañoMatriz):
            for j in range(tamañoMatriz):
                matriz[i][j]=request.POST[[i][j]]
    print(matriz)

    return render(request, './metodosPage/luDirecta.html')

def biseccionView(request):
    datos = ()
    if request.method == 'POST':
        fx = request.POST["funcion"]
        tol = request.POST["tolerancia"]
        Tol = float(tol)
        niter = request.POST["iteraciones"]
        Niter = int(niter)
        xs = request.POST["xs"]
        Xs = float(xs)
        xi = request.POST["xi"]
        Xi = float(xi)

        datos = biseccion(fx, Tol, Niter, Xs, Xi)

    if datos:
        return render(request, './metodosPage/biseccion.html', {'data': datos})

    return render(request, './metodosPage/biseccion.html')

def secanteView(request):
    datos = ()
    if request.method == 'POST':
        fx = request.POST["funcion"]
        tol = request.POST["tolerancia"]
        Tol = float(tol)
        niter = request.POST["iteraciones"]
        Niter = int(niter)
        X0 = request.POST["xs"]
        x0 = float(X0)
        X1 = request.POST["xi"]
        x1 = float(X1)

        datos = secante(fx, Tol, Niter, x0, x1)

    if datos:
        return render(request, './metodosPage/secante.html', {'data': datos})

    return render(request, './metodosPage/secante.html')

def puntoFijoView(request):
    datos = ()
    if request.method == 'POST':
        fx = request.POST["funcion-F"]
        gx = request.POST["funcion-G"]

        x0 = request.POST["vInicial"]
        X0 = float(x0)

        tol = request.POST["tolerancia"]
        Tol = float(tol)

        niter = request.POST["iteraciones"]
        Niter = int(niter)

        datos = puntoFijo(X0,Tol,Niter,fx,gx)
    
    if datos: 
        return render(request, './metodosPage/puntoFijo.html', {'data': datos})


    return render(request, './metodosPage/puntoFijo.html')

def newtonView(request):
    datos = ()
    if request.method == 'POST':
        fx = request.POST["funcion"]
        derf = request.POST["funcion-df"]

        x0 = request.POST["vInicial"]
        X0 = float(x0)

        tol = request.POST["tolerancia"]
        Tol = float(tol)

        niter = request.POST["iteraciones"]
        Niter = int(niter)

        datos = newton(X0, Tol, Niter, fx, derf)
    
    if datos: 
        return render(request, './metodosPage/newton.html', {'data': datos})
    
    return render(request, './metodosPage/newton.html')

def reglaFalsaView(request):
    datos = ()
    if request.method == 'POST':
        fx = request.POST["funcion"]

        x0 = request.POST["lowerinterval"]
        X0 = float(x0)

        xi = request.POST["higherinterval"]
        Xi = float(xi)

        tol = request.POST["tolerancia"]
        Tol = float(tol)

        niter = request.POST["iteraciones"]
        Niter = int(niter)

        datos = reglaFalsa(X0, Xi, Niter, Tol, fx)
        
    
    if datos: 
        return render(request, './metodosPage/reglaFalsa.html', {'data': datos})
    
    return render(request, './metodosPage/reglaFalsa.html')

def raicesMultiplesView(request):
    datos = ()
    if request.method == 'POST':
        Fx = request.POST["funcion"]

        X0 = request.POST["x0"]
        X0 = float(X0)

        niter = request.POST["iteraciones"]
        Niter = int(niter)

        Tol = request.POST["tolerancia"]
        Tol = float(Tol)

        datos = raicesMultiples(Fx,X0,Tol,Niter)
        Errors = datos["errors"]

    if datos:
        return render(request, "./metodosPage/raicesMultiples.html", {"data":datos, "errors": Errors})

    return render(request, './metodosPage/raicesMultiples.html')

def jacobiSeidelView(request):
    datos = []
    if request.method == 'POST':
        mA = toMatrix(request.POST["matrizA"])
        Vx0 = toVector(request.POST["vectorX0"])
        Vb = toVector(request.POST["vectorB"])

        iter = request.POST["iteraciones"]
        Niter = int(iter)

        Tol = request.POST["tolerancia"]
        Tol = float(Tol)

        datos = jacobi(mA,Vb,Vx0,Tol, Niter)

        if datos:
            return render(request, "./metodosPage/jacobi.html", {"data":datos})

    return render(request, './metodosPage/jacobi.html')

def gaussSeidelView(request):
    datos = []
    if request.method == 'POST':
        mA = toMatrix(request.POST["matrizA"])
        Vx0 = toVector(request.POST["vectorX0"])
        Vb = toVector(request.POST["vectorB"])

        niter = request.POST["iteraciones"]
        Niter = int(niter)

        Tol = request.POST["tolerancia"]
        Tol = float(Tol)

        datos = gaussSeidel(mA,Vb,Vx0,Tol, Niter)

        if datos:
            return render(request, "./metodosPage/gaussSeidel.html", {"data":datos})

    return render(request, './metodosPage/gaussSeidel.html')

def sorView(request):
    datos = []
    if request.method == 'POST':
        mA = toMatrix(request.POST["matrizA"])
        Vx0 = toVector(request.POST["vectorX0"])
        Vb = toVector(request.POST["vectorB"])
        print(mA)
        print(Vx0)
        print(Vb)
        w = request.POST["wValue"]
        W = float(w)

        niter = request.POST["iteraciones"]
        Niter = int(niter)

        Tol = request.POST["tolerancia"]
        Tol = float(Tol)

        datos = sor(mA,Vb,Vx0,W,Tol, Niter)

        if datos:
            return render(request, "./metodosPage/sor.html", {"data":datos})

    return render(request, './metodosPage/sor.html')

def splineLinealView(request):
    if request.method == 'POST':
        x = request.POST["x"]
        X = x.split(",")
        X = [float(i) for i in X]
        y = request.POST["y"]
        Y = y.split(",")
        Y = [float(i) for i in Y]
        
        output = splineLineal(X,Y)
    
        Coef = ""
        Tracers = ""
        Errors = output["errors"]

        if(len(Errors)==0):    
            Dic = splineOutput(output)
            data = Dic.split("\n")
            Coef = [data[7], data[8], data[9]]
            Tracers = [data[12], data[13], data[14]]

        return render(request, "./metodosPage/spline-lineal.html",{"coef":Coef, "tracers":Tracers ,"errors":Errors})

    return render(request, "./metodosPage/spline-lineal.html")

def splineCuadraticaView(request):
    if request.method == 'POST':
        x = request.POST["x"]
        X = x.split(",")
        X = [float(i) for i in X]
        y = request.POST["y"]
        Y = y.split(",")
        Y = [float(i) for i in Y]
        
        output = splineCuadratica(X,Y)
    
        Coef = ""
        Tracers = ""
        Errors = output["errors"]

        if(len(Errors)==0):    
            Dic = splineOutput(output)
            data = Dic.split("\n")
            Coef = [data[7], data[8], data[9]]
            Tracers = [data[12], data[13], data[14]]

        return render(request, "./metodosPage/spline-cuadratica.html",{"coef":Coef, "tracers":Tracers ,"errors":Errors})

    return render(request, "./metodosPage/spline-cuadratica.html")

def splineCubicaView(request):
    if request.method == 'POST':
        x = request.POST["x"]
        X = x.split(",")
        X = [float(i) for i in X]
        y = request.POST["y"]
        Y = y.split(",")
        Y = [float(i) for i in Y]
        
        output = splineCubica(X,Y)
    
        Coef = ""
        Tracers = ""
        Errors = output["errors"]

        if(len(Errors)==0):    
            Dic = splineOutput(output)
            data = Dic.split("\n")
            Coef = [data[7], data[8], data[9]]
            Tracers = [data[12], data[13], data[14]]

        return render(request, "./metodosPage/spline-cubica.html",{"coef":Coef, "tracers":Tracers ,"errors":Errors})

    return render(request, "./metodosPage/spline-cubica.html")

def vandermondeView(request):
    if request.method=='POST':
        vectorX=toVector(request.POST["vectorX"])
        vectorY=toVector(request.POST["vectorY"])
        datos=vandermonde(vectorX, vectorY)
        if datos:
            return render(request, "./metodosPage/vandermonde.html", {"data":datos})
    return render(request, "./metodosPage/vandermonde.html")

#Corregir
def newtonDifDevView(request):
    datos = ()
    if request.method == 'POST':
        x = request.POST["x"]
        X = x.split(",")
        X = [float(i) for i in X]
        y = request.POST["y"]
        Y = y.split(",")
        Y = [float(i) for i in Y]

        datos = output = newtonDivDif(X,Y)

    if datos:    
        return render(request, "./metodosPage/newton-DifDev.html",{"datos":datos})

    return render(request, "./metodosPage/newton-DifDev.html")

#Metodos auxiliar
def toMatrix(matrixStr):
    matrixStr = matrixStr.replace(" ","")
    matrixStr = matrixStr.replace("\n","")
    rows = matrixStr.split(";")
    auxM = []
    for row in rows:
        splitedRow = row.split(",")
        auxR = []
        for num in splitedRow:
            auxR.append(float(num))
        auxM.append(auxR)
    return auxM

def toVector(vectorStr):

    splitedVector = vectorStr.split(",")
    auxV = list()
    for num in splitedVector:
        auxV.append(float(num))
    return auxV

def splineOutput(output):
    stringOutput = f'\n"Metodo"\n'
    stringOutput += "\nResults:\n"
    stringOutput += "\nTracer coefficients:\n\n"
    rel = output["results"]
    i = 0
    aux = rel.shape
    while i < aux[0] :
        j = 0
        while j < aux[1]:
            stringOutput += '{:^6f}'.format(rel[i,j]) +"  "
            j += 1
        i += 1
        stringOutput += "\n"
    stringOutput += "\n Tracers:\n"
    i = 0
    while i < aux[0] :
        j = 0
        if aux[1] == 2:
            stringOutput += format(rel[i,0],"6f") +"x"
            stringOutput += format(rel[i,1],"+.6f") 
        elif aux[1] == 3:
            stringOutput += format(rel[i,0],"6f") +"x^2"
            stringOutput += format(rel[i,1],"+.6f") +"x"
            stringOutput += format(rel[i,2],"+.6f")
        elif aux[1] == 4:
            stringOutput += format(rel[i,0],"6f") +"x^3"
            stringOutput += format(rel[i,1],"+.6f") +"x^2"
            stringOutput += format(rel[i,2],"+.6f") +"x"
            stringOutput += format(rel[i,3],"+.6f")

        i += 1
        stringOutput += "\n"
    stringOutput += "\n______________________________________________________________\n"

    return stringOutput

def newtonDifDevOutput(output):
    stringOutput = f'\n"Metodo"\n'
    stringOutput += "\nResults:\n"
    stringOutput += "\nDivided differences table:\n\n"
    rel = output["D"]
    stringOutput += '{:^7f}'.format(rel[0,0]) +"   //L \n"

    stringOutput += "\nNewton's polynomials coefficents:\n\n"
    rel = output["Coef"]
    
    stringOutput += "\nNewton interpolating polynomials:\n\n"
    rel = output["Coef"]
    i = 0
    while i < len(rel) :
        stringOutput += '{:^7f}'.format(rel[i,0]) +"x^3"
        stringOutput += format(rel[i,1],"+.6f") + "x^2"
        stringOutput += format(rel[i,2],"+.6f") + "x" 
        stringOutput += format(rel[i,3],"+.6f") + "   //L \n"
        i += 1

    stringOutput += "\n______________________________________________________________\n"
    return stringOutput