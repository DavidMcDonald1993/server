from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, 
        name='index'),
    path("favicon.ico", views.favicon, 
        name="favicon"),
    path('pass_app', views.index, 
        name='index'),
    path("pass_app/upload/", views.upload_file,
        name="upload"),
    path("pass_app/success/", views.success, 
        name="success",),
]
