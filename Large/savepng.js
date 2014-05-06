
d3.select("#save").on("click", function(){
  var html = d3.select(".intgraph svg")
        .attr("version", 1.1)
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .node().parentNode.innerHTML;
  console.log(html);
  var canvas = document.getElementById("canvas");
  canvg(canvas, html);
  var a = document.createElement("a");
  a.download = "sample.svg";
  a.href = canvas.toDataURL("image/svg");
  a.click();
});
