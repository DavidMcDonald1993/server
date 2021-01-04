function showDiv(button){

    const categories = ["ALL", "PROPERTIES", "ANTITARGETS", "EFFECTS", 
        "GENE_EXPRESSION", "MECHANISMS", "METABOLISM", 
        "TOXICITY", "TRANSPORTERS"];
    
    categories.forEach( category => {

        var style = button.value == category ? 'block' : 'none';

        document.getElementById("nav-" + category).style.cssText = "display:" + style + ";";

    } );


}