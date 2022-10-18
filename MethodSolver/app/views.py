from django.shortcuts import render

# Create your views here.

def menu(request):
    return render(request, 'home.html')

def biseccion(request):
    return render(request, './metodos/biseccion.html')

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

