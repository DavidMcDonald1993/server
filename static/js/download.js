$(document).ready(function() { // dollar sign is jQuery alias
    $('a[data-auto-download]').each(function(){
      var $this = $(this);
      setTimeout(function() {
        window.location = $this.attr('href'); // use href attribute to change page
      }, 10000); // 10 seconds
    });
  });