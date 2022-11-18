from django.http import HttpResponseRedirect
from django.shortcuts import render

from app.methods import *
# Create your views here.


def infoView(request):
    return render(request, 'home.html')


def menuView(request):
    return render(request, 'metodosMenu.html')


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


def reglaFalsa(request):
    return render(request, './metodosPage/reglaFalsa.html')
