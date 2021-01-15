"""npaiengine URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.1/topics/http/urls/
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
from django.urls import path, include

from django.views.generic.base import RedirectView
from django.contrib.staticfiles.storage import staticfiles_storage

from django.views.generic.base import TemplateView

urlpatterns = [
    path("", #include("npaiengine.urls"))
        include("activity_prediction.urls", )),
    path('admin/', admin.site.urls),
    path("activity_prediction/", 
        include("activity_prediction.urls")),
    path("natural_products/", 
        include("natural_products.urls")),
    path("hit_optimisation/",
        include("hit_optimisation.urls")),
    
    
    path("favicon.ico", 
        RedirectView.as_view(url=staticfiles_storage.url("favicon.ico")),),
    path("robots.txt",
        TemplateView.as_view(template_name="robots.txt", content_type="text/plain"),),
    
]
