from django.urls import path, re_path

from . import views

urlpatterns = [


    path("compounds/", 
        views.all_compounds_view, name="all_compounds"),
    path("compounds/all/CNP<str:compound_id>",
        views.compound_info_view, name="compound_info"),


    path("targets",
        views.all_targets_view, name="all_targets"),
    path("targets/all/<str:target>",
        views.target_info_view, name="target_info"),
        
    path('targets/screening/select/', 
        views.target_select_view, name='target_select'),
    path('targets/screening/hits/', 
        views.show_target_hits_view, name='target_hits'),
  

    path("pathways",
        views.all_pathways_view, name="all_pathways"),
    # path("pathways/<str:pathway>:_:<str:organism>",
    path('pathways/screening/select/', 
        views.pathway_select_view, name='pathway_select'),
    path('pathways/screening/hits/', 
        views.show_pathway_hits_view, name='pathway_hits'),
    re_path(r"pathways/all/(?P<pathway_organism>.*)$",
        views.pathway_info_view, name="pathway_info"),


    path("reactions",
        views.all_reactions_view, name="all_reactions"),
    # re_path(r"reactions/(?P<reaction>\d+):_:<str:organism>$",
    path('reactions/screening/select/', 
        views.reaction_select_view, name='reaction_select'),
    path('reactions/screening/hits/', 
        views.show_reaction_hits_view, name='reaction_hits'),

    re_path(r"reactions/all/(?P<reaction_organism>.*)$",
        views.reaction_info_view, name="reaction_info"),

    path('export/', 
        views.export_hits_view, name='export_hits'),
    path("smiles/",
        views.export_hit_smiles_view, name="export_smiles"),
    path('optimise/', 
        views.optimise_target_hits_view, name='optimise_target_hits'),


]