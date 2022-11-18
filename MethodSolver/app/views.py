from django.http import HttpResponseRedirect
from django.shortcuts import render

from app.methods import biseccion
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
        niter = request.POST["iteraciones"]
        xs = request.POST["xs"]
        xi = request.POST["xi"]

        datos = biseccion(fx, tol, niter, xs, xi)
        print(datos)

    if datos:
        print("dentro del if", datos)
        return render(request, './metodosPage/biseccion.html', {'data': datos})
    

    return render(request, './metodosPage/biseccion.html')
    


def gausspiv(request):
    return render(request, 'home.html')


def newton(request):
    return render(request, 'home.html')


def pf(request):
    return render(request, 'home.html')


def rf(request):
    return render(request, 'home.html')


def secante(request):
    return render(request, 'home.html')


def taylorcos(request):
    return render(request, 'home.html')


def taylorsen(request):
    return render(request, 'home.html')
