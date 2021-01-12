from django.urls import path

from . import views

urlpatterns = [
    path('', views.login_page_view, 
        name='login_page'),
    path("logout", views.logout_page_view,
        name="logout"),
    path('login_unsuccessful', views.login_unsuccessful_view, 
        name='login_unsuccessful'),
    path('pass_app', views.index_view, 
        name='index'),
    path("pass_app/upload/", views.upload_file_view,
        name="upload"),
    path("pass_app/success/", views.success_view, 
        name="success",),
        
    path("download/<str:token>", views.download_view, 
        name="download",),
    path("download_error", views.download_error_view, 
        name="download_error",),
]
