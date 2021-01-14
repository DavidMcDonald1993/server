function showDiv(select){

    const categories = ["PROPERTIES", "ALL", "ANTITARGETS", "EFFECTS", 
        "GENE_EXPRESSION", "MECHANISMS", "METABOLISM", 
        "TOXICITY", "TRANSPORTERS", "PATHWAYS", "REACTIONS"];
    
    categories.forEach( category => {

        var style = select.value == category ? 'block' : 'none';

        document.getElementById("nav-" + category).style.cssText = "display:" + style + ";";

    } );

};
