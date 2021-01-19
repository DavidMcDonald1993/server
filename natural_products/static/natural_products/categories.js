currentCategory = "#PROPERTIES";


function showDiv(anchor){

    // const categories = ["PROPERTIES", "ALL", "ANTITARGETS", "EFFECTS", 
    //     "GENE_EXPRESSION", "MECHANISMS", "METABOLISM", 
    //     "TOXICITY", "TRANSPORTERS", "PATHWAYS", "REACTIONS"];

    category = anchor.getAttribute('href');
    // alert(category)
    if (currentCategory != null) {
        document.getElementById("nav-" + currentCategory).style.cssText = "display: none;";
    }
    document.getElementById("nav-" + category).style.cssText = "display: inline;";
    currentCategory = category;
    
    // categories.forEach( category => {

    //     var style = select.value == category ? 'block' : 'none';

    //     document.getElementById("nav-" + category).style.cssText = "display:" + style + ";";

    // } );

};


// $('#categorySidebar a').click(function(e) {
//     // e.preventDefault();

//     category = this.getAttribute('href');
//     alert(category);
// });