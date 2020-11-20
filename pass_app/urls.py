from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, 
        name='index'),
    path("upload/", views.upload_file,
        name="upload"),
    path("success/", views.success, 
        name="success",),
    # path('<str:compound_name>/', 
    #     views.compound_detail, 
    #     name='compound_detail'),
]
