// $(document).ready(makeDataTable(4));

function makeDataTable(table_ids, sort_cols) {

  table_ids.forEach((table_id, idx) => {
    var sort_col = sort_cols[idx]
    $("#" + table_id).dataTable({
      "pageLength": 10,
      "autoWidth": false,
      "fixedColumns": true,
      "order": [sort_col, "desc" ],
    });
  });
 
}