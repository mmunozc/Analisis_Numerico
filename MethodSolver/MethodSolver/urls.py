"""MethodSolver URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from app import views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.infoView, name="infoView"),
    path('menu', views.menuView, name="menuView"),
    path('iterativeMethods', views.iterativeMethodsView, name="iterativeMethodsView"),
    path('splineMethods', views.splineMethodsView, name="splineMethodsView"),
    path('grafica', views.graficaView, name="graficaView"),
    path('biseccion', views.biseccionView, name="biseccion"),
    path('secante', views.secanteView, name="secante"),
    path('puntoFijo', views.puntoFijoView, name="puntoFijo"),
    path('newton', views.newtonView, name="newton"),
    path('reglaFalsa', views.reglaFalsaView, name="reglaFalsa"),
    path('luDirecta', views.luDirectaView, name='luDirecta'),
    path('raicesMultiples', views.raicesMultiplesView, name='raicesMultiples'),
    path('jacobi', views.jacobiSeidelView, name='jacobi'),
    path('gaussSeidel', views.gaussSeidelView, name='gaussSeidel'),
    path('sor', views.sorView, name='sor'),
    path('splineLineal', views.splineLinealView, name='splineLineal'),
    path('splineCuadratica', views.splineCuadraticaView, name='aplineCuadratica'),
    path('splineCubica', views.splineCubicaView, name='splineCubica'),
    path('vandermonde', views.vandermondeView, name='vandermonde'),
    path('newton-difDiv', views.newtonDifDevView, name='newton-difDiv'),
]
