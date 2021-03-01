currentCategory = "#PROPERTIES";

function showDiv(anchor){

    category = anchor.getAttribute('href');
    if (currentCategory != null) {
        document.getElementById("nav-" + currentCategory).style.cssText = "display: none;";
    }
    document.getElementById("nav-" + category).style.cssText = "display: inline;";
    currentCategory = category;
    
};