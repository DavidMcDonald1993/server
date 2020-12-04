from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('targets/', views.target, name='target'),
    path('results/', 
        views.results, name='results'),
    path("compounds/", 
        views.all_compounds, name="all_compounds"),
    path("compounds/CNP<str:compound_id>",
        views.compound_info, name="compound_info"),
    path('download/', 
        views.download, name='download'),
    path('optimise/', 
        views.optimise, name='optimise'),
]