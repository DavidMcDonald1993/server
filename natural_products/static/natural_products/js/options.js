currentCategory = null;

function showDiv(option){

    category = option.value.replace(" ", "_");
    console.log(category);
    if (currentCategory != null) {
        document.getElementById("div-" + currentCategory).style.cssText = "display: none;";
    }
    document.getElementById("div-" + category).style.cssText = "display: inline;";
    currentCategory = category;
    
};