from django.urls import path, re_path

from .views import compound_views, target_views, pathway_views, reaction_views, uniprot_views, export

urlpatterns = [

    path("compounds/", 
        compound_views.all_compounds_view, name="all_compounds"),
    path("compounds/all/CNP<str:compound_id>",
        compound_views.compound_info_view, name="compound_info"),

    path("targets",
        target_views.all_targets_view, name="all_targets"),
    path("targets/all/<str:target>",
        target_views.target_info_view, name="target_info"),
    path('targets/screening/select/', 
        target_views.target_select_view, name='target_select'),
    path('targets/screening/hits/', 
        target_views.show_target_hits_view, name='target_hits'),
  
    path("pathways",
        pathway_views.all_pathways_view, name="all_pathways"),
    path('pathways/screening/select/', 
        pathway_views.pathway_select_view, name='pathway_select'),
    path('pathways/screening/hits/', 
        pathway_views.show_pathway_hits_view, name='pathway_hits'),
    re_path(r"pathways/all/(?P<pathway_organism>.*)$",
        pathway_views.pathway_info_view, name="pathway_info"),

    path("reactions",
        reaction_views.all_reactions_view, name="all_reactions"),
    path('reactions/screening/select/', 
        reaction_views.reaction_select_view, name='reaction_select'),
    path('reactions/screening/hits/', 
        reaction_views.show_reaction_hits_view, name='reaction_hits'),
    re_path(r"reactions/all/(?P<reaction_organism>.*)$",
        reaction_views.reaction_info_view, name="reaction_info"),
    
    path("uniprot",
        uniprot_views.all_accs_view, name="all_uniprot"),
    path("uniprot/all/<str:acc>",
        uniprot_views.acc_info_view, name="reaction_info"),
    # re_path(r"reactions/(?P<reaction>\d+):_:<str:organism>$",
    # path('reactions/screening/select/', 
    #     views.reaction_select_view, name='reaction_select'),
    # path('reactions/screening/hits/', 
    #     views.show_reaction_hits_view, name='reaction_hits'),


    path('export/', 
        export.export_hits_view, name='export_hits'),
    path("smiles/",
        export.export_hit_smiles_view, name="export_smiles"),
    path('optimise/', 
        export.optimise_target_hits_view, name='optimise_target_hits'),


]