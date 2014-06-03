// menu list for visualization

var result = [];


(function(){
  // toggle for motion
  var moact = function(){

    $(".motion").show();
  }
  var modeact = function(){
    $(".motion").hide();
  }


  // toggle for highlight protein targets
  var targetact = function(){
    $(".stat").show();
  }
  var targetdeact = function(){
    $(".stat").hide();
  }

  var menuitems = ["Disable Movement", "Highlight Target Proteins", "Brush Statistics"];
  var funcselect = [ moact, targetact ];
  var funcunselect = [ modeact, targetdeact ];

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
        $(".menu li").each(function(){
          $(this).addClass('ui-selected');
        });
      },
      selected: function( e, u ) {
        var idx = menuitems.indexOf(u.selected.innerText);
        funcselect[idx]();
      },
      unselected: function( e, u ) {
        var idx = menuitems.indexOf(u.unselected.innerText);
        funcunselect[idx]();
      },
      stop: function() {
        $( ".ui-selected", this ).each(function() {
          selected_name = this.innerText;
          result.push( selected_name );
        });
      }
    });
  });
})()
