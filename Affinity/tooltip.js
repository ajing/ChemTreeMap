// Tool tip for each node

var maketooltip = function(edata){
    console.log(edata);
    var image_src = "../Image/" + edata.name;
    var html = '<input type="button" value="As Reference" onclick="saverefer(\'' + image_src + '\');">';
    html = html + "<h3>" + edata.name + "</h3>";
    console.log(edata.name);
    html = html + "<ul>";
    for (var key in edata.continuous){
      html = html + "<li>" + key + ": " + edata.continuous[key] + "</li>";
    }
    for (var key in edata.nominal){
      html = html + "<li>" + key + ": " + edata.nominal[key] + "</li>";
    }
    html = html + "</ul>";
    html = html + '<img id="molimage" src="' + image_src + '">';
    return html;
};

window.saverefer = function(img_src){
    d3.select("#molimage").attr("src", img_src);
};

var tooltipconfig =
{
  items: ".node",
  content: function(){
    var eledata = d3.select(this).data()[0];
    return maketooltip(eledata);
  },
  position: { my: "left top-15" },
  show: {
    effect: "slideDown",
    delay: 50,
  },
  open: function(event, ui)
  {
      if (typeof(event.originalEvent) === 'undefined')
      {
          return false;
      }

      var $id = $(ui.tooltip).attr('id');

      // close any lingering tooltips
      $('div.ui-tooltip').not('#' + $id).remove();

      // ajax function to pull in data and add it to the tooltip goes here
  },
  close: function(event, ui)
  {
      ui.tooltip.hover(function()
      {
          $(this).stop(true).fadeTo(400, 1);
      },
      function()
      {
          $(this).fadeOut('400', function()
          {
              $(this).remove();
          });
      });
  }
};

var moltip = $("svg").tooltip(tooltipconfig);
