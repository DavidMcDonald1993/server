function makeDataTable(args) {

  if (args.length == 2) {
    table_ids = args[0]
    sort_cols = args[1]
    table_ids.forEach((table_id, idx) => {
      var sort_col = sort_cols[idx]
      $("#" + table_id).dataTable({
        "pageLength": 10,
        "autoWidth": false,
        "fixedColumns": true,
        "order": [],
      });
    });
  } else {
    $("table.display").dataTable({
      "pageLength": 10,
      "autoWidth": false,
      "fixedColumns": true,
      "order": [],
    });

  }

 
}