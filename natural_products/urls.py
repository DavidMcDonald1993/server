from django.urls import path

from . import views

urlpatterns = [
    path('', views.index_view, name='index'),
 
    path('target_selection/', 
        views.target_select_view, name='target_select'),
    path('target_hits/', 
        views.show_target_hits_view, name='target_hits'),
    path('download/', 
        views.download_hits_view, name='download_hits'),
    path("smiles/",
        views.download_hit_smiles_view, name="download_smiles"),
    path('optimise/', 
        views.optimise_target_hits_view, name='optimise_target_hits'),
 
    path('pathway_selection/', 
        views.pathway_select_view, name='pathway_select'),
    path('pathway_hits/', 
        views.show_pathway_hits_view, name='pathway_hits'),

    path('reaction_selection/', 
        views.reaction_select_view, name='reaction_select'),
    path('reaction_hits/', 
        views.show_reaction_hits_view, name='reaction_hits'),

    path("compounds/", 
        views.all_compounds_view, name="all_compounds"),
    path("compounds/CNP<str:compound_id>",
        views.compound_info_view, name="compound_info"),

    # path("pathway_enrichment", 
    #     views.pathway_enrichment, name="enrichment")

]