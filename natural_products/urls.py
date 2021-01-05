from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
 
    path('target_selection/', 
        views.target_select, name='target_select'),
    path('target_hits/', 
        views.show_target_hits, name='target_hits'),
    path('download/', 
        views.download_target_hits, name='download_target_hits'),
    path('optimise/', 
        views.optimise_target_hits, name='optimise_target_hits'),
 
    path('pathway_selection/', 
        views.pathway_select, name='pathway_select'),
    path('pathway_hits/', 
        views.show_pathway_hits, name='pathway_hits'),

    path('reaction_selection/', 
        views.reaction_select, name='reaction_select'),
    path('reaction_hits/', 
        views.show_reaction_hits, name='reaction_hits'),

    path("compounds/", 
        views.all_compounds, name="all_compounds"),
    path("compounds/CNP<str:compound_id>",
        views.compound_info, name="compound_info"),

    # path("pathway_enrichment", 
    #     views.pathway_enrichment, name="enrichment")

]