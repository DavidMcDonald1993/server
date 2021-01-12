from django.urls import path

from hit_optimisation import views

urlpatterns = [
    path('', views.index_view, 
        name='index'),
    path("upload/", views.upload_view, 
        name='upload'),
    path("success/", views.success_view, 
        name="success")
]
