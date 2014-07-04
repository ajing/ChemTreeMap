// Tool tip for each node

var maketooltip = function(edata){
    var html = "<h3>" + edata.name + "</h3>";
    html = html + "<ul>";
    for (var key in edata.continuous){
      html = html + "<li>" + key + ": " + edata.continuous[key] + "</li>";
    }
    for (var key in edata.nominal){
      html = html + "<li>" + key + ": " + edata.nominal[key] + "</li>";
    }
    html = html + "</ul>";
    var image_src = "../Image/" + edata.name;
    html = html + '<img id="molimage" src="' + image_src + '">';
    return html;
}

$(document).tooltip({
  items:".node",
  content: function(){
    var eledata = d3.select(this).data()[0];
    return maketooltip(eledata);
  },
  show: {
    effect: "slideDown",
    delay: 400,
  }
})
