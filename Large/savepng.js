
d3.select("#save").on("click", function(){
  var html = d3.select(".intgraph svg")
        .attr("version", 1.1)
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .node().parentNode.innerHTML;
  //console.log(html);
  var imgsrc = 'data:image/svg+xml;base64,'+ btoa(html);
  var canvas = document.querySelector("canvas"),
  context = canvas.getContext("2d");
  var image = new Image;
  image.src = imgsrc;
  image.onload = function() {
  context.drawImage(image, 0, 0);
  //save and serve it as an actual filename
  binaryblob();
  var a = document.createElement("a");
  a.download = "sample.png";
  a.href = canvas.toDataURL("image/png");
  a.click();
  };
});

function binaryblob(){
    var byteString = atob(document.querySelector("canvas").toDataURL().replace(/^data:image\/(png|jpg);base64,/, "")); //wtf is atob?? https://developer.mozilla.org/en-US/docs/Web/API/Window.atob
    var ab = new ArrayBuffer(byteString.length);
    var ia = new Uint8Array(ab);
    for (var i = 0; i < byteString.length; i++) {
        ia[i] = byteString.charCodeAt(i);
    }
    var dataView = new DataView(ab);
    var blob = new Blob([dataView], {type: "image/png"});
    var DOMURL = self.URL || self.webkitURL || self;
    var newurl = DOMURL.createObjectURL(blob);
    var img = '<img src="'+newurl+'">';
    d3.select("#img").html(img);
}

