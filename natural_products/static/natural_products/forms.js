$(document).ready(function() {
    $('select').change(function() {

      // alert("press");

      var empty = false;
      $('select').each(function() {
          if ($(this).val() == '') {
              empty = true;
          }
      });
      if (empty) {
          $('#submit').attr('disabled', 'disabled'); 
      } else {
          $('#submit').removeAttr('disabled');
      }
    });
});