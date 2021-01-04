from django.urls import path

from . import views

urlpatterns = [
    path('', views.login_page, 
        name='login_page'),
    path('login_unsuccessful', views.login_unsuccessful, 
        name='login_unsuccessful'),
    path("favicon.ico", views.favicon, 
        name="favicon"),
    path('pass_app', views.index, 
        name='index'),
    path("pass_app/upload/", views.upload_file,
        name="upload"),
    path("pass_app/success/", views.success, 
        name="success",),
]
