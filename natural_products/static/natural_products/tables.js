
$(document).ready(function(){
  $('table.display.table.data').removeAttr("width").dataTable({
    "pageLength": 10,
    "autoWidth": false,
    "fixedColumns": true
  });
});

function setDataTablesVisible() {
  var tables = document.getElementsByClassName("table display table data")

  console.log(tables.length)

  // tables.forEach(table => {
  for (var i=0; i<tables.length; i++){
    table = tables[i]
    table.style.cssText = "display:inline"
  };
}