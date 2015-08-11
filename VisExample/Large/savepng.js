$(function() {
  $( "button[id=savepng]" )
    .button()
    .click(function( event ) {
        var html = d3.select(".intgraph svg")
              .attr("version", 1.1)
              .attr("xmlns", "http://www.w3.org/2000/svg")
              .node().parentNode.innerHTML;
        $.post("downloadpng", {svg : html}, function(){
          var a = document.createElement("a");
          a.download = "sample.png";
          a.href = "../sample.png";
          a.click();
        });
        //console.log(html);
        //var canvas = document.getElementById("canvas");
        //canvg(canvas, html);
        //var a = document.createElement("a");
        //a.download = "sample.png";
        //a.href = canvas.toDataURL("image/png");
        //a.click();
    });
});

$(function() {
  $( "button[id=savesvg]" )
    .button()
    .click(function( event ) {
        var html = d3.select(".intgraph svg")
              .attr("version", 1.1)
              .attr("xmlns", "http://www.w3.org/2000/svg")
              .node().parentNode.innerHTML;
        $.post("downloadsvg", {svg : html}, function(){
          var a = document.createElement("a");
          a.download = "sample.svg";
          a.href = "../sample.svg";
          a.click();
        });
    });
});
