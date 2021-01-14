function showHidden() {
    var hiddenElements = document.getElementsByClassName("hiddenDiv");

    for (var i=0; i<hiddenElements.length; i++){ // forEach not valid
        hiddenElements[i].style.cssText = "display:inline"
    }
}
  
$(document).ready(showHidden);
  