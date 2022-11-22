from django.http import HttpResponseRedirect
from django.shortcuts import render
import numpy as np
from app.methods import *
# Create your views here.


def infoView(request):
    return render(request, 'home.html')


def menuView(request):
    return render(request, 'metodosMenu.html')

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

#def luDirectaView(request):



def puntoFijoView(request):
    datos = ()
    if request.method == 'POST':
        fx = request.POST["funcion-F"]
        gx = request.POST["funcion-f-prima"]

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

def raphsonView(request):
    datos = ()
    if request.method == 'POST':
        f = request.POST["funcion-F"]
        fp = request.POST["funcion-F-prima"]

        xi = request.POST["vInicial"]
        Xi = float(xi)

        tol = request.POST["tolerancia"]
        Tol = float(tol)

        niter = request.POST["iteraciones"]
        Niter = int(niter)



        datos = raphson(f, fp, Xi, Tol ,Niter)
    
    if datos: 
        return render(request, './metodosPage/raphson.html', {'data': datos})


    return render(request, './metodosPage/raphson.html')