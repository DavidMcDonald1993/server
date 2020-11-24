
$(document).ready(function(){
  $('#myTable').dataTable(
    {"pageLength": 25}
  );
})

function setTableVisible() {
  document.getElementById("myTable").style.cssText = "width:100%;display:block;"
}