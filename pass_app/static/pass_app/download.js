$(function() {
    $('a[data-auto-download]').each(function(){
        var $this = $(this);
        setTimeout(function() {
            window.location = $this.attr('href');
         }, 1000);
    });
});