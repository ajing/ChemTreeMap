// menu list for visualization

(function(){
  // toggle for motion
  var moact = function(){
    node.call(drag);
    node.on('click', click);
    $(".motion").show();
  }
  var modeact = function(){
    node.on('mousedown.drag', null);
    node.on('click', null);
    $(".motion").hide();
  }

  // toggle for highlight protein targets
  var targetact = function(){
    $(".target").show();
  }
  var targetdeact = function(){
    $(".target").hide();
  }

  // toggle for brushing histogram
  var brushact = function(){
    $(".histogram").show();
    graph.on('mousedown.zoom', null);
    brushing();
  }
  var brushdeact = function(){
    $(".histogram").hide();
    d3.selectAll(".brush").remove();
    node.classed("selected", false);
    graph.call(zoomer);
  }

  var menuitems = ["Enable Movement", "Highlight Target Proteins", "Brush Statistics"];
  var funcselect = [ moact, targetact, brushact ];
  var funcunselect = [ modeact, targetdeact, brushdeact ];

  var table = d3.select('.menu').append('ul')
      .selectAll('li')
      .data(menuitems)
  .enter()
      .append('li')
      .text(function(d) { return d});

  $(function() {
    $( ".menu" ).bind("mousedown", function(e){
      e.metaKey = true;
    }).selectable({
      filter:'li',
      create: function() {
        $(".menu li").each(function(e){
          if (menuitems[e].indexOf("Brush") == -1)
             $(this).addClass('ui-selected');
        });
      },
      selected: function( e, u ) {
        var idx = menuitems.indexOf(u.selected.innerHTML);
        funcselect[idx]();
      },
      unselected: function( e, u ) {
        var idx = menuitems.indexOf(u.unselected.innerHTML);
        funcunselect[idx]();
      }
    });
  });
})()
