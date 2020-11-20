from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('targets/', views.target, name='target'),
    # path('results/?target=<str:target>&threshold=<str:threshold>', 
    path('results/', 
        views.results, name='results'),
    path("compounds/CNP<str:compound_id>",
        views.compound_info, name="compound_info")
]