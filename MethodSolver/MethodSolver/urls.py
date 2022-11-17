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
    path('biseccion', views.biseccionView, name="biseccion"),
    #path('gausspiv', views.menu, name="gausspiv"),
    #path('newton', views.menu, name="newton"),
    #path('pf', views.menu, name="pf"),
    #path('rf', views.menu, name="rf"),
    #path('secante', views.menu, name="secante"),
    #path('taylorcos', views.menu, name="taylorcos"),
    #path('taylorsen', views.menu, name="taylorsen"),

]
