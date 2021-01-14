function showDiv(select){

    const categories = ["ANTITARGETS", "EFFECTS", 
        "GENE_EXPRESSION", "MECHANISMS", "METABOLISM", 
        "TOXICITY", "TRANSPORTERS"];
    
    categories.forEach( category => {

        var style = select.value == category ? 'block' : 'none';

        document.getElementById("div-" + category).style.cssText = "display:" + style + ";";

    } );

};
