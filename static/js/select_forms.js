$(document).ready(function() {
    $('select').change(function() {

      var empty = false;
      $('select').each(function() {
          var element = $(this);
          if (element.css("display") != "none" && element.val() == '') {
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