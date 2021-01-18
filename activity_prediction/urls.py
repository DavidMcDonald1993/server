from django.urls import path

from . import views

urlpatterns = [
    path('', views.index_view, 
        name='index'),
    path('activity_prediction', views.index_view, 
        name='index'),

    path('login', views.login_view, 
        name='login_page'),
    path("logout", views.logout_view,
        name="logout"),
    path('login_unsuccessful', views.login_unsuccessful_view, 
        name='login_unsuccessful'),

    path("activity_prediction/upload/", views.upload_file_view,
        name="upload"),
    path("activity_prediction/success/", views.success_view, 
        name="success",),
        
    path("download/<str:token>", views.download_view, 
        name="download",),
    path("download_error", views.download_error_view, 
        name="download_error",),

]
